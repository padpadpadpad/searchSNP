//! Checking whether the search sequences also occur elsewhere in the reference genome.
//!
//! A search sequence is usually unique to the focal window, but not always: if the region
//! is duplicated in the reference (or a repetitive element matches it), reads from the
//! other copy(s) are counted too, which silently inflates the ref/alt proportions.
//!
//! This module indexes the reference once per window length and, for every search
//! sequence, records
//!
//!  * `n_occ_ref`: how many times it occurs in the whole reference (any chromosome,
//!    either strand — the same canonical k-mer semantics as the read counter),
//!  * `elsewhere`: the positions of the occurrences OUTSIDE the focal window.
//!
//! Occurrences are matched with the canonical (either-strand) form and reported as the
//! forward-strand start position of the matching window.
//!
//! The index is a flat sorted array of (canonical k-mer, packed position) pairs, about
//! 16 bytes per base of genome, built in one pass (a few hundred ms for a 5 Mb genome).

use crate::kmers::{base_code, canonical_of};
use crate::reference::Reference;
use crate::sequences::VariantSeqs;
use crate::vcf::Variant;
use std::collections::HashMap;

/// Maximum number of elsewhere-positions kept per sequence (the count stays complete;
/// only the printed list is capped, by the writers as well).
pub const MAX_KEPT_POSITIONS: usize = 1000;

/// One (canonical k-mer, position) entry of the index, sorted by k-mer then position.
struct Entry {
    key: u64,
    /// chromosome index in the high 32 bits, 1-based forward position in the low 32 bits
    occ: u64,
}

#[inline]
fn pack(chrom: usize, pos: usize) -> u64 {
    ((chrom as u64) << 32) | pos as u64
}

#[inline]
fn unpack(occ: u64) -> (usize, usize) {
    ((occ >> 32) as usize, (occ & 0xffff_ffff) as usize)
}

/// All (canonical k-mer, position) pairs of one window length, sorted.
struct OccIndex {
    k: usize,
    entries: Vec<Entry>,
}

impl OccIndex {
    fn build(reference: &Reference, k: usize) -> OccIndex {
        assert!((1..=31).contains(&k), "k must be between 1 and 31");
        let mut entries: Vec<Entry> = Vec::new();
        let mask: u64 = (1u64 << (2 * k)) - 1;
        let rc_shift = 2 * (k - 1);
        for (ci, c) in reference.chroms.iter().enumerate() {
            // rolling canonical k-mer, the same arithmetic as the read counter
            let (mut fwd, mut rc, mut valid) = (0u64, 0u64, 0usize);
            for (i, &b) in c.seq.iter().enumerate() {
                match base_code(b) {
                    Some(code) => {
                        fwd = ((fwd << 2) | code) & mask;
                        rc = (rc >> 2) | ((3 - code) << rc_shift);
                        valid += 1;
                        if valid >= k {
                            // 1-based position of the window's first base
                            entries.push(Entry { key: fwd.min(rc), occ: pack(ci, i + 2 - k) });
                        }
                    }
                    // N or other symbol breaks the window
                    None => valid = 0,
                }
            }
        }
        entries.sort_unstable_by_key(|e| (e.key, e.occ));
        OccIndex { k, entries }
    }

    /// Every position (chrom index, 1-based) whose canonical k-mer equals that of `seq`,
    /// sorted by chromosome then position.
    fn occurrences(&self, seq: &[u8]) -> Vec<(usize, usize)> {
        if seq.len() != self.k {
            return Vec::new();
        }
        let Some(key) = canonical_of(seq) else { return Vec::new() };
        let lo = self.entries.partition_point(|e| e.key < key);
        let hi = self.entries.partition_point(|e| e.key <= key);
        self.entries[lo..hi].iter().map(|e| unpack(e.occ)).collect()
    }
}

/// Annotate every search sequence with its reference occurrences (see the module docs).
///
/// Fills `n_occ_ref` / `elsewhere` on each `SearchSeq`; occurrences at the sequence's own
/// focal window count towards `n_occ_ref` but not towards `elsewhere`. Returns the number
/// of sequences that occur somewhere else in the reference.
///
/// One index per distinct window length is built on first use (there is almost always
/// only one length; more appear when short chromosomes force truncated windows).
pub fn annotate(seqs: &mut [VariantSeqs], vars: &[Variant], reference: &Reference) -> usize {
    let mut cache: HashMap<usize, OccIndex> = HashMap::new();
    let mut n_elsewhere = 0;
    for (list, v) in seqs.iter_mut().zip(vars) {
        for s in &mut list.seqs {
            let Some(_) = canonical_of(&s.seq) else { continue }; // non-ACGT: not searchable
            let idx = cache.entry(s.seq.len()).or_insert_with(|| OccIndex::build(reference, s.seq.len()));
            let occ = idx.occurrences(&s.seq);
            s.n_occ_ref = occ.len() as u32;
            s.elsewhere = occ
                .into_iter()
                .filter(|&(c, p)| c != v.chrom || p != s.window_start)
                .take(MAX_KEPT_POSITIONS)
                .collect();
            if !s.elsewhere.is_empty() {
                n_elsewhere += 1;
            }
        }
    }
    n_elsewhere
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::Chrom;
    use crate::sequences::build_all;
    use crate::vcf::Variant;

    /// A genome with only one 7-mer duplicated: "ACGTTGC" at 1:1 and again at 1:21.
    /// (Chosen so no other 7-mer, nor the reverse complement of the ones used in the
    /// tests, appears more than once.)
    const REF2: &str = "ACGTTGCAAGCTTACGGATCACGTTGCTCAGGTCA";

    fn var(pos: usize, r: &str, a: &str) -> Variant {
        Variant { chrom: 0, pos, r: r.as_bytes().to_vec(), a: a.as_bytes().to_vec(), ref_matches: true }
    }

    fn annotate_one(seq: &[u8], window_start: usize) -> (u32, Vec<(usize, usize)>) {
        let mut vs = VariantSeqs {
            seqs: vec![crate::sequences::SearchSeq {
                seq: seq.to_vec(),
                is_ref: true,
                window_start,
                n_occ_ref: 0,
                elsewhere: Vec::new(),
            }],
            window_bp: seq.len(),
            status: crate::sequences::SearchStatus::Full,
        };
        let vars = [var(4, "T", "A")];
        annotate(std::slice::from_mut(&mut vs), &vars, &ref2());
        let s = &vs.seqs[0];
        (s.n_occ_ref, s.elsewhere.clone())
    }

    fn ref2() -> Reference {
        Reference::from_chroms(vec![Chrom { name: "x".into(), seq: REF2.as_bytes().to_vec() }])
    }

    #[test]
    fn duplicated_sequence_is_flagged_with_positions() {
        // "ACGTTGC" occurs at 1:1 and 1:21; window_start = 1, so one is elsewhere
        let (n, elsewhere) = annotate_one(b"ACGTTGC", 1);
        assert_eq!(n, 2);
        assert_eq!(elsewhere, vec![(0, 21)]);
    }

    #[test]
    fn unique_sequence_has_no_elsewhere() {
        // "AAGCTTA" occurs once, at 1:8 — the focal window here
        let (n, elsewhere) = annotate_one(b"AAGCTTA", 8);
        assert_eq!(n, 1);
        assert!(elsewhere.is_empty());
    }

    #[test]
    fn absent_sequence_counts_zero() {
        let (n, elsewhere) = annotate_one(b"ACACACA", 1);
        assert_eq!((n, elsewhere.is_empty()), (0, true));
    }

    #[test]
    fn reverse_strand_occurrence_is_found() {
        // The genome contains "AAGCTTA" forward at 1:8 (but not its revcomp). Because
        // matching is canonical (either strand), searching the revcomp "TAAGCTT" still
        // reports the forward occurrence at 1:8.
        let mut rc = b"AAGCTTA".to_vec();
        rc.reverse();
        for b in rc.iter_mut() {
            *b = match *b { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => *b };
        }
        assert_eq!(rc, b"TAAGCTT");
        assert!(!REF2.contains("TAAGCTT")); // absent as a forward match
        let (n, elsewhere) = annotate_one(&rc, 1);
        assert_eq!(n, 1);
        assert_eq!(elsewhere, vec![(0, 8)]);
    }

    #[test]
    fn end_to_end_through_build_all() {
        // variant at 1:4 (flank 3): ref window starts at 1 and equals "ACGTTGC",
        // which also occurs at 1:21 -> ref sequences flagged, alt sequences absent.
        let r = ref2();
        let vars = [var(4, "T", "A")];
        let mut all = build_all(&vars, &r, 3, 1 << 20, 1);
        let n = annotate(&mut all, &vars, &r);
        assert_eq!(n, 1); // only the ref-labelled sequence occurs elsewhere
        let refseq = all[0].seqs.iter().find(|s| s.is_ref).unwrap();
        assert_eq!(refseq.n_occ_ref, 2);
        assert_eq!(refseq.elsewhere, vec![(0, 21)]);
        let altseq = all[0].seqs.iter().find(|s| !s.is_ref).unwrap();
        assert_eq!(altseq.n_occ_ref, 0);
        assert!(altseq.elsewhere.is_empty());
    }
}
