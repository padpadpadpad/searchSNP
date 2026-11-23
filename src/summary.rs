//! Turning search sequences + k-mer counts into the per-VCF-row summary table.

use crate::kmers::{canonical_of, KmerIndex};
use crate::reference::Reference;
use crate::sequences::VariantSeqs;
use crate::vcf::Variant;
use std::collections::HashMap;
use std::io::Write;

/// For one variant: the k-mer index entries (distinct) of its ref- and alt-labelled
/// haplotypes, plus bookkeeping.
pub struct VariantKmers {
    /// distinct k-mer indices of haplotypes where the focal variant is REF
    pub ref_idx: Vec<usize>,
    /// distinct k-mer indices of haplotypes where the focal variant is ALT
    pub alt_idx: Vec<usize>,
    /// k-mers that occur under BOTH labels of this variant (e.g. a long insertion fills the
    /// window so a nearby focal deletion cannot change the string). Reads of these cannot
    /// tell ref from alt, so they are excluded from `ref_idx`, `alt_idx` and the summary.
    pub ambiguous: Vec<usize>,
    /// number of search sequences that could be searched (made only of A/C/G/T)
    pub n_search_seqs: usize,
    /// true if a k-mer is also used by a variant at another position
    pub shared: bool,
}

/// Collect all search sequences into one `KmerIndex` and record, for each variant,
/// which k-mers belong to it. Sequences containing N or other non-ACGT symbols cannot
/// be searched and are left out.
pub fn index_sequences(seqs: &[VariantSeqs], vars: &[Variant]) -> (KmerIndex, Vec<VariantKmers>) {
    let mut index = KmerIndex::new();
    // first (chrom, pos) that used each k-mer; a different later user means "shared"
    let mut owner: HashMap<usize, (usize, usize)> = HashMap::new();
    let mut out: Vec<VariantKmers> = Vec::with_capacity(seqs.len());

    for (v, list) in vars.iter().zip(seqs) {
        let mut ref_idx: Vec<usize> = Vec::new();
        let mut alt_idx: Vec<usize> = Vec::new();
        let mut n = 0;
        let mut shared = false;
        for s in &list.seqs {
            let Some(c) = canonical_of(&s.seq) else { continue };
            n += 1;
            let i = index.insert(c, s.seq.len());
            let target = if s.is_ref { &mut ref_idx } else { &mut alt_idx };
            if !target.contains(&i) {
                target.push(i);
            }
            // used by a variant at a different position?
            let site = (v.chrom, v.pos);
            if *owner.entry(i).or_insert(site) != site {
                shared = true;
            }
        }
        // a k-mer under both labels would count the same reads as ref and alt: set it aside
        let ambiguous: Vec<usize> = ref_idx.iter().copied().filter(|i| alt_idx.contains(i)).collect();
        ref_idx.retain(|i| !ambiguous.contains(i));
        alt_idx.retain(|i| !ambiguous.contains(i));
        out.push(VariantKmers { ref_idx, alt_idx, ambiguous, n_search_seqs: n, shared });
    }
    // a k-mer first claimed by an earlier variant at another position flags only the later
    // one above; flag the earlier owners too so the flag is symmetric
    let mut users: HashMap<usize, Vec<usize>> = HashMap::new();
    for (vi, vk) in out.iter().enumerate() {
        for &i in vk.ref_idx.iter().chain(vk.alt_idx.iter()).chain(vk.ambiguous.iter()) {
            users.entry(i).or_default().push(vi);
        }
    }
    for list in users.values() {
        let first = (vars[list[0]].chrom, vars[list[0]].pos);
        if list.iter().any(|&vi| (vars[vi].chrom, vars[vi].pos) != first) {
            for &vi in list {
                out[vi].shared = true;
            }
        }
    }
    (index, out)
}

/// Sum of the counts of the given k-mer indices.
fn sum_counts(idx: &[usize], counts: &[u64]) -> u64 {
    idx.iter().map(|&i| counts[i]).sum()
}

/// Write the summary table: one row per VCF row (per ALT allele), in VCF order.
///
/// `fwd` / `rev` are the k-mer counts from the forward / reverse reads. Missing k-mers count as 0, and
/// `prop_alt` is 0 when there are no reads at all (`total` = 0).
///
/// `ref_check`: whether the reference-occurrence check ran (see `refcheck`). If not, the
/// `shared_in_reference` column is left empty instead of false.
pub fn write_summary<W: Write>(
    mut w: W,
    vars: &[Variant],
    kmers: &[VariantKmers],
    seqs: &[VariantSeqs],
    reference: &Reference,
    fwd: &[u64],
    rev: &[u64],
    ref_check: bool,
) -> std::io::Result<()> {
    writeln!(
        w,
        "chrom,pos,ref,alt,\
         ref_count_fwd,alt_count_fwd,total_fwd,prop_alt_fwd,\
         ref_count_rev,alt_count_rev,total_rev,prop_alt_rev,\
         window_bp,ref_matches,multiallelic,n_alts_at_site,\
         n_search_seqs,search_status,shared_kmers,n_ambiguous_seqs,shared_in_reference"
    )?;

    // number of rows (ALT alleles) at each chrom/pos, to flag multi-allelic sites
    let mut per_site: HashMap<(usize, usize), usize> = HashMap::new();
    for v in vars {
        *per_site.entry((v.chrom, v.pos)).or_insert(0) += 1;
    }

    for ((v, k), sq) in vars.iter().zip(kmers).zip(seqs) {
        let n_alts = per_site[&(v.chrom, v.pos)];
        // shared_in_reference: a searched sequence also occurs outside its focal window
        let shared_ref = if ref_check {
            sq.seqs.iter().any(|s| !s.elsewhere.is_empty()).to_string()
        } else {
            String::new()
        };
        write!(
            w,
            "{},{},{},{}",
            reference.chroms[v.chrom].name,
            v.pos,
            String::from_utf8_lossy(&v.r),
            String::from_utf8_lossy(&v.a)
        )?;
        // the same four columns for forward and reverse reads
        for counts in [fwd, rev] {
            let r = sum_counts(&k.ref_idx, counts);
            let a = sum_counts(&k.alt_idx, counts);
            let total = r + a;
            let prop = if total == 0 { 0.0 } else { a as f64 / total as f64 };
            write!(w, ",{r},{a},{total},{prop:.6}")?;
        }
        // diagnostics; n_ambiguous = sequences set aside because they occur under both labels
        writeln!(
            w,
            ",{},{},{},{},{},{},{},{},{}",
            sq.window_bp,
            v.ref_matches,
            n_alts > 1,
            n_alts,
            k.n_search_seqs,
            sq.status.as_str(),
            k.shared,
            k.ambiguous.len(),
            shared_ref
        )?;
    }
    Ok(())
}

/// Write the per-sequence table: one row per searched sequence, in VCF order, with the
/// read counts of that sequence (either strand) in the forward (`count_fwd`) and reverse (`count_rev`) reads.
///
/// A variant that was skipped by one of our rules (see `search_status`) has no sequences;
/// it gets one placeholder row with an empty sequence and empty counts, so that every VCF
/// row appears. A searched sequence that was never seen has counts of 0.
///
/// `ref_check`: whether the reference-occurrence check ran (see `refcheck`). If not, the
/// three `ref_occurrences*` columns are left empty instead of 0.
pub fn write_sequence_counts<W: Write>(
    mut w: W,
    vars: &[Variant],
    seqs: &[VariantSeqs],
    kmers: &[VariantKmers],
    index: &KmerIndex,
    reference: &Reference,
    fwd: &[u64],
    rev: &[u64],
    ref_check: bool,
) -> std::io::Result<()> {
    writeln!(
        w,
        "chrom,pos,ref,alt,type,sequence,count_fwd,count_rev,window_bp,search_status,shared_kmers,ambiguous,\
         ref_occurrences,ref_occurrences_elsewhere,ref_elsewhere_positions"
    )?;
    for ((v, vs), vk) in vars.iter().zip(seqs).zip(kmers) {
        // the first columns are identical for every row of this variant
        let head = format!(
            "{},{},{},{}",
            reference.chroms[v.chrom].name,
            v.pos,
            String::from_utf8_lossy(&v.r),
            String::from_utf8_lossy(&v.a)
        );
        if vs.seqs.is_empty() {
            // skipped variant: placeholder row, nothing was searched
            writeln!(w, "{head},,,,,{},{},{},,,,", vs.window_bp, vs.status.as_str(), vk.shared)?;
            continue;
        }
        for s in &vs.seqs {
            // counts of this sequence's canonical k-mer; empty if it cannot be searched (N)
            let canon = canonical_of(&s.seq);
            let hit = canon.and_then(|c| index.get(c, s.seq.len()));
            // ambiguous = same k-mer is also under the other label of this variant
            let amb = hit.map_or(false, |i| vk.ambiguous.contains(&i));
            let (cf, cr) = match hit {
                Some(i) => (fwd[i].to_string(), rev[i].to_string()),
                None => (String::new(), String::new()),
            };
            // reference occurrences; empty when not checked or the sequence is not searchable
            let (n_occ, n_else, pos_list) = match (ref_check, canon) {
                (true, Some(_)) => (
                    s.n_occ_ref.to_string(),
                    s.elsewhere.len().to_string(),
                    s.elsewhere
                        .iter()
                        .take(5)
                        .map(|&(c, p)| format!("{}:{}", reference.chroms[c].name, p))
                        .collect::<Vec<_>>()
                        .join(";"),
                ),
                _ => (String::new(), String::new(), String::new()),
            };
            writeln!(
                w,
                "{head},{},{},{cf},{cr},{},{},{},{amb},{n_occ},{n_else},{pos_list}",
                if s.is_ref { "ref" } else { "alt" },
                String::from_utf8_lossy(&s.seq),
                vs.window_bp,
                vs.status.as_str(),
                vk.shared
            )?;
        }
    }
    Ok(())
}
