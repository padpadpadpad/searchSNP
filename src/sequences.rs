//! Building the search sequences for every variant.
//!
//! For each focal variant we build every haplotype sequence of length `2*flank + 1`
//! that could be seen in reads around it:
//!
//!   [`flank` reference bases before POS] + [edited sequence from POS onwards] ... cut to length
//!
//! "Edited" means any combination of the focal allele and the alleles of other variants
//! that start within `flank` bases of the focal position. Each sequence is labelled
//! `ref` or `alt` depending ONLY on whether the focal variant carries its REF or ALT.
//!
//! Indels: VCF indels start with an anchor base at POS. We keep `flank` bases before POS,
//! then the ALT allele (anchor + inserted bases), then reference bases after the REF
//! allele until the sequence is `2*flank + 1` long. For a deletion the reference after the
//! deleted bases is pulled in; for a long insertion only the left junction fits.
//!
//! Chromosome ends: if there are fewer than `flank` bases before POS (or after it), the
//! window is SHIFTED along the chromosome so it keeps its full length (the variant is then
//! off-centre but the sequence is still exact and as specific as any other). Only if the
//! whole chromosome is shorter than the window is a shorter window used ("truncated"), and
//! if even that is shorter than `MIN_WINDOW` the variant is skipped.

use crate::reference::Reference;
use crate::vcf::Variant;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::thread;

/// One haplotype sequence to search for.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SearchSeq {
    /// The sequence (forward strand, upper case), exactly `2*flank + 1` bases.
    pub seq: Vec<u8>,
    /// True if the focal variant carries its REF allele in this haplotype.
    pub is_ref: bool,
    /// 1-based id of the haplotype background this sequence belongs to, within its
    /// variant. The REF and ALT sequences of one background share it, so a ref/alt pair
    /// can be identified without relying on row order.
    pub haplotype: u32,
    /// 1-based reference position the window was built from (its left end, in the
    /// unedited reference coordinate frame; may differ from POS at chromosome ends or
    /// when a deletion makes the window slide).
    pub window_start: usize,
    /// How many times the sequence occurs in the whole reference (any chromosome, either
    /// strand), filled in by `refcheck::annotate`; 0 also means "not checked" or the
    /// sequence is not searchable (contains N).
    pub n_occ_ref: u32,
    /// Reference positions (chromosome index, 1-based) of occurrences OUTSIDE the focal
    /// window, also from `refcheck::annotate`.
    pub elsewhere: Vec<(usize, usize)>,
}

/// Smallest window (bp) that will be searched. For flanks so small that the full window
/// is already shorter than this (flank < 7), the full window length is the minimum.
pub const MIN_WINDOW: usize = 15;

/// What happened to the window of a variant (reported in the `search_status` column).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SearchStatus {
    /// Normal window: `flank` bases either side, full length.
    Full,
    /// Too close to a chromosome end to centre the window; slid along to keep full length.
    Shifted,
    /// Chromosome shorter than the full window; the whole-chromosome-length window is used.
    Truncated,
    /// Window would be shorter than the minimum; nothing searched.
    SkippedTooShort,
    /// More haplotype combinations than `--max-combos`; nothing searched.
    SkippedTooManyCombos,
}

impl SearchStatus {
    /// Text written to the output table.
    pub fn as_str(self) -> &'static str {
        match self {
            SearchStatus::Full => "full",
            SearchStatus::Shifted => "shifted",
            SearchStatus::Truncated => "truncated",
            SearchStatus::SkippedTooShort => "skipped_too_short",
            SearchStatus::SkippedTooManyCombos => "skipped_too_many_combos",
        }
    }
}

/// Everything built for one variant.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VariantSeqs {
    pub seqs: Vec<SearchSeq>,
    /// Length in bp of the sequences searched (0 if skipped).
    pub window_bp: usize,
    pub status: SearchStatus,
}

/// One possible edit at a position: replace `reflen` reference bases with `allele`.
struct Edit<'a> {
    allele: &'a [u8],
    reflen: usize,
}

/// Build the search sequences of variant `i`.
///
/// `by_chrom` holds, per chromosome, variant indices sorted by position (used to find
/// neighbours quickly with a binary search).
fn build_variant(
    i: usize,
    vars: &[Variant],
    by_chrom: &[Vec<usize>],
    reference: &Reference,
    flank: usize,
    max_combos: u128,
) -> VariantSeqs {
    let f = &vars[i];
    let chrom = &reference.chroms[f.chrom];
    let seq = &chrom.seq;
    let want_len = 2 * flank + 1;
    let skipped = |status| VariantSeqs { seqs: Vec::new(), window_bp: 0, status };

    // Window length actually used: the full length unless the chromosome is shorter.
    let eff_len = want_len.min(seq.len());
    if eff_len < MIN_WINDOW.min(want_len) {
        eprintln!("warning: chromosome too short for a window at {}:{}, variant skipped", chrom.name, f.pos);
        return skipped(SearchStatus::SkippedTooShort);
    }

    // Neighbours are variants whose POS is within +/- flank of the focal POS.
    let lo = f.pos.saturating_sub(flank);
    let hi = f.pos + flank;
    let list = &by_chrom[f.chrom];
    let first = list.partition_point(|&k| vars[k].pos < lo);

    // For each position, the list of possible EDITS (non-reference choices).
    // Choice 0 at every position is "leave the reference as it is".
    // The focal position only ever offers the focal ALT: other variants at the
    // focal position are ignored, as they cannot be on the same haplotype.
    let mut sites: BTreeMap<usize, Vec<Edit>> = BTreeMap::new();
    sites.insert(f.pos, vec![Edit { allele: &f.a, reflen: f.r.len() }]);
    for &k in &list[first..] {
        let v = &vars[k];
        if v.pos > hi {
            break;
        }
        if v.pos == f.pos {
            continue;
        }
        let edits = sites.entry(v.pos).or_default();
        // skip exact duplicates (same variant listed twice)
        if !edits.iter().any(|e| e.allele == &v.a[..] && e.reflen == v.r.len()) {
            edits.push(Edit { allele: &v.a, reflen: v.r.len() });
        }
    }

    let positions: Vec<usize> = sites.keys().copied().collect();
    let edit_lists: Vec<&Vec<Edit>> = sites.values().collect();
    let n_sites = positions.len();
    let focal_slot = positions.iter().position(|&p| p == f.pos).unwrap();

    // number of haplotype combinations = product of (edits + 1 for "no edit") per site
    let total: u128 = edit_lists.iter().map(|e| e.len() as u128 + 1).product();
    if total > max_combos {
        eprintln!(
            "warning: {total} haplotype combinations at {}:{} exceeds --max-combos, variant skipped",
            chrom.name, f.pos
        );
        return skipped(SearchStatus::SkippedTooManyCombos);
    }

    // The window normally starts `flank` bases before POS. Too close to the chromosome
    // start (pos <= flank) it starts at base 1 instead, i.e. it is shifted.
    let start0 = lo.max(1);
    let mut shifted = f.pos <= flank;

    // `choice[s]` = 0 for "no edit" at site s, or 1.. for edit number choice-1.
    // It is advanced like an odometer (first site fastest), like R's expand.grid.
    let mut choice = vec![0usize; n_sites];
    let mut out: Vec<SearchSeq> = Vec::new();
    let mut seen: HashSet<(Vec<u8>, bool)> = HashSet::new();
    // background haplotype (the `choice` vector with the focal slot zeroed) -> 1-based id
    let mut haplotypes: HashMap<Vec<usize>, u32> = HashMap::new();
    let mut next_haplotype = 1u32;

    for _ in 0..total {
        // collect the edits chosen for this combination, in position order
        let mut applied: Vec<(usize, &Edit)> = Vec::new();
        for s in 0..n_sites {
            if choice[s] > 0 {
                applied.push((positions[s], &edit_lists[s][choice[s] - 1]));
            }
        }

        // two edits that overlap on the reference cannot both be on one haplotype: skip
        let mut conflict = false;
        let mut prev_end = 0usize; // last reference position covered by the previous edit
        for (p, e) in &applied {
            if *p <= prev_end {
                conflict = true;
                break;
            }
            prev_end = p + e.reflen - 1;
        }

        if !conflict {
            // Build `eff_len` bases from `start`. If the chromosome ends first (a deletion
            // or the right-hand edge), slide the start left by the missing amount and
            // rebuild, until the sequence is full or the start reaches base 1.
            let mut start = start0;
            let buf = loop {
                let mut buf: Vec<u8> = Vec::with_capacity(eff_len + 16);
                let mut cur = start; // next reference position (1-based) to copy
                for (p, e) in &applied {
                    copy_ref(&mut buf, seq, cur, *p);
                    buf.extend_from_slice(e.allele);
                    cur = p + e.reflen;
                    if buf.len() >= eff_len {
                        break;
                    }
                }
                if buf.len() < eff_len {
                    let need = eff_len - buf.len();
                    copy_ref(&mut buf, seq, cur, cur + need);
                }
                buf.truncate(eff_len);
                if buf.len() == eff_len || start == 1 {
                    break buf;
                }
                start = start.saturating_sub(eff_len - buf.len()).max(1);
                shifted = true;
            };

            // a haplotype that cannot reach full length even from base 1 is dropped
            if buf.len() == eff_len {
                let is_ref = choice[focal_slot] == 0;
                if seen.insert((buf.clone(), is_ref)) {
                    // the pair is identified by the other sites' choices: the focal slot
                    // is what distinguishes REF from ALT within one background
                    let mut background = choice.clone();
                    background[focal_slot] = 0;
                    let haplotype = *haplotypes.entry(background).or_insert_with(|| {
                        let id = next_haplotype;
                        next_haplotype += 1;
                        id
                    });
                    out.push(SearchSeq {
                        seq: buf,
                        is_ref,
                        haplotype,
                        window_start: start,
                        n_occ_ref: 0,
                        elsewhere: Vec::new(),
                    });
                }
            }
        }

        // advance the odometer
        for s in 0..n_sites {
            choice[s] += 1;
            if choice[s] <= edit_lists[s].len() {
                break;
            }
            choice[s] = 0;
        }
    }

    let status = if eff_len < want_len {
        SearchStatus::Truncated
    } else if shifted {
        SearchStatus::Shifted
    } else {
        SearchStatus::Full
    };
    VariantSeqs { seqs: out, window_bp: eff_len, status }
}

/// Append reference bases [from, to_excl) (1-based, to_excl exclusive) to `buf`,
/// stopping at the end of the chromosome.
fn copy_ref(buf: &mut Vec<u8>, seq: &[u8], from: usize, to_excl: usize) {
    let a = from - 1;
    let b = (to_excl - 1).min(seq.len());
    if a < b {
        buf.extend_from_slice(&seq[a..b]);
    }
}

/// Build search sequences for all variants, using `threads` threads.
/// The result has one entry per variant, in the same order as `vars`.
pub fn build_all(
    vars: &[Variant],
    reference: &Reference,
    flank: usize,
    max_combos: u128,
    threads: usize,
) -> Vec<VariantSeqs> {
    // per-chromosome variant indices sorted by position
    let mut by_chrom: Vec<Vec<usize>> = vec![Vec::new(); reference.chroms.len()];
    for (i, v) in vars.iter().enumerate() {
        by_chrom[v.chrom].push(i);
    }
    for l in by_chrom.iter_mut() {
        l.sort_by_key(|&k| vars[k].pos);
    }

    let n = vars.len();
    if n == 0 {
        return Vec::new();
    }
    // split the variants into contiguous chunks, one per thread, so order is preserved
    let nt = threads.max(1).min(n);
    let chunk = (n + nt - 1) / nt;
    thread::scope(|s| {
        let handles: Vec<_> = (0..n)
            .step_by(chunk)
            .map(|st| {
                let en = (st + chunk).min(n);
                let by_chrom = &by_chrom;
                s.spawn(move || {
                    (st..en)
                        .map(|i| build_variant(i, vars, by_chrom, reference, flank, max_combos))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles.into_iter().flat_map(|h| h.join().unwrap()).collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::Chrom;

    // Test reference, 1-based positions:
    //  1 A  2 C  3 G  4 T  5 T  6 G  7 C  8 A  9 A 10 G
    // 11 C 12 T 13 T 14 A 15 C 16 G 17 G 18 A 19 T 20 C
    const REF: &str = "ACGTTGCAAGCTTACGGATC";

    fn reference() -> Reference {
        Reference::from_chroms(vec![Chrom { name: "c".into(), seq: REF.as_bytes().to_vec() }])
    }

    fn var(pos: usize, r: &str, a: &str) -> Variant {
        Variant {
            chrom: 0,
            pos,
            r: r.as_bytes().to_vec(),
            a: a.as_bytes().to_vec(),
            ref_matches: true,
        }
    }

    /// Build for variant `focal` among `vars`, and return (ref seqs, alt seqs) as sorted strings.
    fn run(vars: &[Variant], focal: usize, flank: usize) -> (Vec<String>, Vec<String>) {
        let r = reference();
        let all = build_all(vars, &r, flank, 1 << 20, 1);
        let mut refs = Vec::new();
        let mut alts = Vec::new();
        for s in &all[focal].seqs {
            let t = String::from_utf8(s.seq.clone()).unwrap();
            if s.is_ref { refs.push(t) } else { alts.push(t) }
        }
        refs.sort();
        alts.sort();
        (refs, alts)
    }

    #[test]
    fn lone_snp() {
        // window 5..11 = TGC [A] AGC
        let (r, a) = run(&[var(8, "A", "G")], 0, 3);
        assert_eq!(r, vec!["TGCAAGC"]);
        assert_eq!(a, vec!["TGCGAGC"]);
    }

    #[test]
    fn deletion_pulls_in_reference_after() {
        // AAG (8-10) -> A: TGC + A + reference from 11 (CTT) = TGCACTT
        let (r, a) = run(&[var(8, "AAG", "A")], 0, 3);
        assert_eq!(r, vec!["TGCAAGC"]);
        assert_eq!(a, vec!["TGCACTT"]);
    }

    #[test]
    fn short_insertion() {
        // A -> ATT: TGC + ATT + reference from 9 (AGC...) cut to 7 = TGCATTA
        let (_, a) = run(&[var(8, "A", "ATT")], 0, 3);
        assert_eq!(a, vec!["TGCATTA"]);
    }

    #[test]
    fn long_insertion_keeps_left_junction_only() {
        // insert is longer than the flank: only part of it fits, length stays 7
        let (_, a) = run(&[var(8, "A", "ACCCCCCC")], 0, 3);
        assert_eq!(a, vec!["TGCACCC"]);
    }

    #[test]
    fn two_snps_give_four_haplotypes_labelled_by_focal() {
        let vars = [var(8, "A", "G"), var(10, "G", "T")];
        let (r, a) = run(&vars, 0, 3);
        assert_eq!(r, vec!["TGCAAGC", "TGCAATC"]);
        assert_eq!(a, vec!["TGCGAGC", "TGCGATC"]);
    }

    #[test]
    fn haplotype_ids_pair_the_ref_and_alt_of_one_background() {
        // lone SNP: one background, so both labels share id 1
        let all = build_all(&[var(8, "A", "G")], &reference(), 3, 1 << 20, 1);
        let ids: Vec<(u32, bool)> = all[0].seqs.iter().map(|s| (s.haplotype, s.is_ref)).collect();
        assert_eq!(ids, vec![(1, true), (1, false)]);

        // two sites: two backgrounds, each with one ref and one alt sequence
        let vars = [var(8, "A", "G"), var(10, "G", "T")];
        let all = build_all(&vars, &reference(), 3, 1 << 20, 1);
        let ids: Vec<(u32, bool)> = all[0].seqs.iter().map(|s| (s.haplotype, s.is_ref)).collect();
        assert_eq!(ids, vec![(1, true), (1, false), (2, true), (2, false)]);
    }

    #[test]
    fn snp_with_upstream_deletion_keeps_full_length() {
        // GC (6-7) -> G upstream shifts the focal base but the length stays 7
        let vars = [var(8, "A", "G"), var(6, "GC", "G")];
        let (r, a) = run(&vars, 0, 3);
        assert_eq!(r, vec!["TGAAGCT", "TGCAAGC"]);
        assert_eq!(a, vec!["TGCGAGC", "TGGAGCT"]);
    }

    #[test]
    fn overlapping_alts_never_combined() {
        // deletions 8-10 and 9-10 overlap: they cannot be on the same haplotype
        let vars = [var(8, "AAG", "A"), var(9, "AG", "A")];
        let (r, a) = run(&vars, 0, 3);
        // focal ref: plain reference, or the other deletion on its own
        assert_eq!(r, vec!["TGCAACT", "TGCAAGC"]);
        // focal alt: only the focal deletion (the overlapping one conflicts)
        assert_eq!(a, vec!["TGCACTT"]);
    }

    #[test]
    fn other_variant_at_same_position_is_ignored() {
        let vars = [var(8, "A", "G"), var(8, "A", "C")];
        let (r, a) = run(&vars, 0, 3);
        assert_eq!(r, vec!["TGCAAGC"]);
        assert_eq!(a, vec!["TGCGAGC"]);
    }

    #[test]
    fn duplicate_rows_do_not_duplicate_sequences() {
        let vars = [var(8, "A", "G"), var(8, "A", "G")];
        let (r, a) = run(&vars, 0, 3);
        assert_eq!((r.len(), a.len()), (1, 1));
    }

    /// Status and window length of variant `focal` on a given reference sequence.
    fn status_on(chrom: &str, vars: &[Variant], focal: usize, flank: usize) -> (SearchStatus, usize, Vec<String>) {
        let r = Reference::from_chroms(vec![Chrom { name: "c".into(), seq: chrom.as_bytes().to_vec() }]);
        let all = build_all(vars, &r, flank, 1 << 20, 1);
        let seqs = all[focal].seqs.iter().map(|s| String::from_utf8(s.seq.clone()).unwrap()).collect();
        (all[focal].status, all[focal].window_bp, seqs)
    }

    #[test]
    fn interior_variant_is_full() {
        let (st, bp, _) = status_on(REF, &[var(8, "A", "G")], 0, 3);
        assert_eq!((st, bp), (SearchStatus::Full, 7));
        // first/last positions that still have a full centred window
        assert_eq!(status_on(REF, &[var(4, "T", "A")], 0, 3).0, SearchStatus::Full);
        assert_eq!(status_on(REF, &[var(17, "G", "A")], 0, 3).0, SearchStatus::Full);
    }

    #[test]
    fn variant_near_start_is_shifted_and_keeps_full_length() {
        // pos 2, flank 3: window slides to start at base 1 instead of base -1
        let (st, bp, seqs) = status_on(REF, &[var(2, "C", "T")], 0, 3);
        assert_eq!((st, bp), (SearchStatus::Shifted, 7));
        assert_eq!(seqs, vec!["ACGTTGC", "ATGTTGC"]);
        // pos 1
        let (st, _, seqs) = status_on(REF, &[var(1, "A", "T")], 0, 3);
        assert_eq!(st, SearchStatus::Shifted);
        assert_eq!(seqs, vec!["ACGTTGC", "TCGTTGC"]);
    }

    #[test]
    fn variant_near_end_is_shifted_and_keeps_full_length() {
        // pos 19, flank 3: window ends at the last base (20) and starts at 14
        let (st, bp, seqs) = status_on(REF, &[var(19, "T", "A")], 0, 3);
        assert_eq!((st, bp), (SearchStatus::Shifted, 7));
        assert_eq!(seqs, vec!["ACGGATC", "ACGGAAC"]);
        // last base
        let (st, _, seqs) = status_on(REF, &[var(20, "C", "G")], 0, 3);
        assert_eq!(st, SearchStatus::Shifted);
        assert_eq!(seqs, vec!["ACGGATC", "ACGGATG"]);
    }

    #[test]
    fn deletion_near_end_still_full_length() {
        // TC (19-20) -> T at the very end: the alt is one base shorter, so the window
        // slides left to stay 7 long
        let (_, _, seqs) = status_on(REF, &[var(19, "TC", "T")], 0, 3);
        assert_eq!(seqs, vec!["ACGGATC", "TACGGAT"]);
    }

    #[test]
    fn chromosome_shorter_than_window_is_truncated() {
        // 20 bp chromosome, default-size window of 23: use the whole chromosome (20 >= 15)
        let (st, bp, seqs) = status_on(REF, &[var(10, "G", "T")], 0, 11);
        assert_eq!((st, bp), (SearchStatus::Truncated, 20));
        assert_eq!(seqs, vec!["ACGTTGCAAGCTTACGGATC", "ACGTTGCAATCTTACGGATC"]);
    }

    #[test]
    fn chromosome_below_minimum_is_skipped() {
        // 10 bp chromosome, flank 11 -> window would be 10 < 15
        let (st, bp, seqs) = status_on("ACGTTGCAAG", &[var(5, "T", "A")], 0, 11);
        assert_eq!((st, bp), (SearchStatus::SkippedTooShort, 0));
        assert!(seqs.is_empty());
        // 15 bp chromosome is just enough
        let (st, bp, _) = status_on("ACGTTGCAAGCTTAC", &[var(5, "T", "A")], 0, 11);
        assert_eq!((st, bp), (SearchStatus::Truncated, 15));
    }

    #[test]
    fn length_is_always_two_flank_plus_one() {
        for flank in [1usize, 5, 8] {
            let (r, a) = run(&[var(10, "G", "T"), var(12, "TT", "T")], 0, flank);
            for s in r.iter().chain(a.iter()) {
                assert_eq!(s.len(), 2 * flank + 1);
            }
        }
    }

    #[test]
    fn dense_cluster_is_capped() {
        let vars = [var(8, "A", "G"), var(9, "A", "C"), var(10, "G", "T")];
        let r = reference();
        // 2*2*2 = 8 combinations; cap at 4 skips the variant
        let all = build_all(&vars, &r, 3, 4, 1);
        assert!(all[0].seqs.is_empty());
        assert_eq!(all[0].status, SearchStatus::SkippedTooManyCombos);
    }

    #[test]
    fn threads_give_same_result_as_one_thread() {
        let vars = [var(8, "A", "G"), var(10, "G", "T"), var(12, "T", "C")];
        let r = reference();
        assert_eq!(build_all(&vars, &r, 3, 1 << 20, 1), build_all(&vars, &r, 3, 1 << 20, 3));
    }
}
