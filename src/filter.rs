//! Removing ref/alt haplotype pairs whose ALT sequence is present in the ancestral reads.
//!
//! A pair is one haplotype background: the REF and ALT sequences that share the same
//! combination of neighbouring alleles (same `SearchSeq::haplotype`). If the ALT member is
//! found in the ancestral reads, the variant was already there before the experiment, so it
//! says nothing about a new mutation. Both members of the pair are removed, and the
//! variant's counts are rebuilt from the pairs that remain.
//!
//! Removing the REF member too is deliberate: keeping it would divide the novel ALT reads
//! by REF reads from a background where the focal ALT was already segregating. A pair is
//! never half-removed.

use crate::kmers::{canonical_of, KmerIndex};
use crate::sequences::{SearchSeq, VariantSeqs};
use crate::summary::{collect_kmers, VariantKmers};
use crate::vcf::Variant;
use std::collections::{HashMap, HashSet};

/// What the ancestral filter did to one VCF row (the `ancestor_status` column).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RowFilter {
    /// No pair was removed.
    Kept,
    /// Some pairs were removed; the counts are rebuilt from the rest.
    Partial,
    /// Every pair was removed, so nothing is counted (the measurement columns are blank).
    AllDropped,
    /// The variant was skipped before filtering, so there were no pairs to assess.
    NotAssessed,
}

impl RowFilter {
    /// Text written to the output table (empty for `NotAssessed`).
    pub fn as_str(self) -> &'static str {
        match self {
            RowFilter::Kept => "kept",
            RowFilter::Partial => "partial",
            RowFilter::AllDropped => "all_dropped",
            RowFilter::NotAssessed => "",
        }
    }
}

/// Counts of pairs kept and removed, for the console summary.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct FilterStats {
    pub pairs: usize,
    pub dropped: usize,
    pub partial: usize,
    pub all_dropped: usize,
}

/// Sum of a sequence's ancestral counts (both directions); 0 if it cannot be searched.
fn ancestor_total(s: &SearchSeq, index: &KmerIndex, fwd: &[u64], rev: &[u64]) -> u64 {
    match canonical_of(&s.seq).and_then(|c| index.get(c, s.seq.len())) {
        Some(i) => fwd[i] + rev[i],
        None => 0,
    }
}

/// Remove every haplotype pair whose ALT sequence occurs at least `min_alt_count` times in
/// the ancestral reads, and rebuild the per-variant k-mer lists from what remains.
///
/// `seqs`, `vars` and `index` are the unfiltered run's. The returned k-mer lists refer to
/// the same index entries, so the existing read counts can be used unchanged.
pub fn filter_pairs(
    seqs: &[VariantSeqs],
    vars: &[Variant],
    index: &KmerIndex,
    ancestor_fwd: &[u64],
    ancestor_rev: &[u64],
    min_alt_count: u64,
) -> (Vec<VariantSeqs>, Vec<VariantKmers>, Vec<RowFilter>, FilterStats) {
    let mut kept_seqs: Vec<VariantSeqs> = Vec::with_capacity(seqs.len());
    let mut row_filter: Vec<RowFilter> = Vec::with_capacity(seqs.len());
    let mut stats = FilterStats::default();

    for vs in seqs {
        // per haplotype: how often was its ALT member seen in the ancestor?
        // (deduplication guarantees at most one ALT sequence per haplotype)
        let mut alt_total: HashMap<u32, u64> = HashMap::new();
        for s in &vs.seqs {
            if !s.is_ref {
                alt_total.insert(s.haplotype, ancestor_total(s, index, ancestor_fwd, ancestor_rev));
            }
        }
        let dropped: HashSet<u32> = alt_total
            .iter()
            .filter(|(_, total)| **total >= min_alt_count)
            .map(|(&hap, _)| hap)
            .collect();

        stats.pairs += alt_total.len();
        stats.dropped += dropped.len();

        // keep the original order, which is stable and comparable with the unfiltered file
        let kept: Vec<SearchSeq> = vs
            .seqs
            .iter()
            .filter(|s| !dropped.contains(&s.haplotype))
            .cloned()
            .collect();

        let status = if vs.seqs.is_empty() {
            RowFilter::NotAssessed
        } else if kept.is_empty() {
            stats.all_dropped += 1;
            RowFilter::AllDropped
        } else if kept.len() < vs.seqs.len() {
            stats.partial += 1;
            RowFilter::Partial
        } else {
            RowFilter::Kept
        };

        kept_seqs.push(VariantSeqs { seqs: kept, window_bp: vs.window_bp, status: vs.status });
        row_filter.push(status);
    }

    let kmers = collect_kmers(index, &kept_seqs, vars);
    (kept_seqs, kmers, row_filter, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequences::SearchStatus;
    use crate::summary::index_sequences;

    fn var(pos: usize, r: &str, a: &str) -> Variant {
        Variant { chrom: 0, pos, r: r.as_bytes().to_vec(), a: a.as_bytes().to_vec(), ref_matches: true }
    }

    fn seq(s: &str, is_ref: bool, haplotype: u32) -> SearchSeq {
        SearchSeq {
            seq: s.as_bytes().to_vec(),
            is_ref,
            haplotype,
            window_start: 1,
            n_occ_ref: 0,
            elsewhere: Vec::new(),
        }
    }

    fn vs(seqs: Vec<SearchSeq>) -> VariantSeqs {
        VariantSeqs { seqs, window_bp: 5, status: SearchStatus::Full }
    }

    /// Counts vector with the given sequence set to `n` (all other sequences 0).
    fn counts_with(index: &KmerIndex, seq: &str, n: u64) -> Vec<u64> {
        let mut c = vec![0u64; index.len()];
        let i = index.get(canonical_of(seq.as_bytes()).unwrap(), seq.len()).unwrap();
        c[i] = n;
        c
    }

    #[test]
    fn drops_the_pair_whose_alt_is_in_the_ancestor() {
        let vars = [var(10, "A", "G")];
        let seqs = [vs(vec![
            seq("ACGTA", true, 1),
            seq("ACGTC", false, 1),
            seq("ACGTT", true, 2),
            seq("ACGAG", false, 2),
        ])];
        let (index, _) = index_sequences(&seqs, &vars);
        let anc_fwd = counts_with(&index, "ACGAG", 3);
        let anc_rev = vec![0u64; index.len()];

        let (kept, kmers, status, stats) =
            filter_pairs(&seqs, &vars, &index, &anc_fwd, &anc_rev, 1);

        assert_eq!(status, [RowFilter::Partial]);
        assert_eq!((stats.pairs, stats.dropped, stats.partial), (2, 1, 1));
        // only the first background survives, both of its members
        let seqs_kept: Vec<(&[u8], bool, u32)> =
            kept[0].seqs.iter().map(|s| (&s.seq[..], s.is_ref, s.haplotype)).collect();
        assert_eq!(seqs_kept, [(b"ACGTA".as_slice(), true, 1), (b"ACGTC".as_slice(), false, 1)]);
        // the k-mer lists point at the same index entries as before
        let idx = |s: &str| index.get(canonical_of(s.as_bytes()).unwrap(), s.len()).unwrap();
        assert_eq!(kmers[0].ref_idx, [idx("ACGTA")]);
        assert_eq!(kmers[0].alt_idx, [idx("ACGTC")]);
        assert_eq!(kmers[0].n_search_seqs, 2);
    }

    #[test]
    fn threshold_keeps_pairs_below_the_minimum() {
        let vars = [var(10, "A", "G")];
        let seqs = [vs(vec![seq("ACGTA", true, 1), seq("ACGTC", false, 1)])];
        let (index, _) = index_sequences(&seqs, &vars);
        let anc_fwd = counts_with(&index, "ACGTC", 1);
        let anc_rev = vec![0u64; index.len()];

        let (kept, _, status, _) = filter_pairs(&seqs, &vars, &index, &anc_fwd, &anc_rev, 2);
        assert_eq!(status, [RowFilter::Kept]);
        assert_eq!(kept[0].seqs.len(), 2);
    }

    #[test]
    fn all_pairs_dropped_leaves_an_empty_row() {
        let vars = [var(10, "A", "G")];
        let seqs = [vs(vec![seq("ACGTA", true, 1), seq("ACGTC", false, 1)])];
        let (index, _) = index_sequences(&seqs, &vars);
        let anc_fwd = counts_with(&index, "ACGTC", 1);
        let anc_rev = vec![0u64; index.len()];

        let (kept, kmers, status, stats) =
            filter_pairs(&seqs, &vars, &index, &anc_fwd, &anc_rev, 1);
        assert_eq!(status, [RowFilter::AllDropped]);
        assert_eq!((stats.pairs, stats.dropped, stats.all_dropped), (1, 1, 1));
        assert!(kept[0].seqs.is_empty());
        assert_eq!((kept[0].window_bp, kept[0].status), (5, SearchStatus::Full));
        assert!(kmers[0].ref_idx.is_empty() && kmers[0].alt_idx.is_empty());
        assert_eq!(kmers[0].n_search_seqs, 0);
    }

    #[test]
    fn unpaired_sequences_are_judged_on_their_own() {
        let vars = [var(10, "A", "G")];
        // a lone REF (its ALT was lost to a conflict) and a lone ALT (its REF was a duplicate)
        let seqs = [vs(vec![seq("ACGTAA", true, 1), seq("ACGTAC", false, 2)])];
        let (index, _) = index_sequences(&seqs, &vars);
        let anc_fwd = counts_with(&index, "ACGTAC", 1);
        let anc_rev = vec![0u64; index.len()];

        let (kept, _, status, _) = filter_pairs(&seqs, &vars, &index, &anc_fwd, &anc_rev, 1);
        assert_eq!(status, [RowFilter::Partial]);
        // the lone ALT is dropped, the lone REF is kept
        assert_eq!(kept[0].seqs.len(), 1);
        assert!(kept[0].seqs[0].is_ref);
        assert_eq!(kept[0].seqs[0].seq, b"ACGTAA");
    }

    #[test]
    fn skipped_variants_are_not_assessed() {
        let vars = [var(10, "A", "G")];
        let seqs = [vs(Vec::new())];
        let (index, _) = index_sequences(&seqs, &vars);
        let anc_fwd = vec![0u64; index.len()];
        let anc_rev = vec![0u64; index.len()];

        let (kept, _, status, stats) = filter_pairs(&seqs, &vars, &index, &anc_fwd, &anc_rev, 1);
        assert_eq!(status, [RowFilter::NotAssessed]);
        assert_eq!(stats, FilterStats::default());
        assert!(kept[0].seqs.is_empty());
    }
}
