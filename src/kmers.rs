//! Exact k-mer counting in FASTQ reads.
//!
//! Replaces `jellyfish count -C ... --if search.fasta reads.fastq.gz`:
//!  * k-mers are stored 2 bits per base in a `u64` (so k <= 31 here);
//!  * k-mers are CANONICAL: a k-mer and its reverse complement are the same key
//!    (the numerically smaller of the two), so reads from either strand match;
//!  * only k-mers in the search set are counted; everything else is ignored;
//!  * a k-mer must match EXACTLY. Reads containing N (or any non-ACGT) break the window.

use flate2::read::MultiGzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::hash::{BuildHasherDefault, Hasher};
use std::io::{BufRead, BufReader};
use std::sync::mpsc::sync_channel;
use std::sync::{Arc, Mutex};
use std::thread;

/// A very fast hasher for u64 keys (multiplicative hashing). The default SipHash is
/// needlessly slow here: keys are already well-mixed integers and we are not
/// exposed to hostile input.
#[derive(Default)]
pub struct FastHasher(u64);

impl Hasher for FastHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        // generic fallback (not used for u64 keys)
        for &b in bytes {
            self.0 = (self.0.rotate_left(5) ^ b as u64).wrapping_mul(0x517c_c1b7_2722_0a95);
        }
    }
    fn write_u64(&mut self, n: u64) {
        // multiply by a large odd constant, then fold the high bits down
        let h = n.wrapping_mul(0x9E37_79B9_7F4A_7C15);
        self.0 = h ^ (h >> 32);
    }
}

type FastMap = HashMap<u64, usize, BuildHasherDefault<FastHasher>>;

/// 2-bit code of a base (A=0, C=1, G=2, T=3), or None for anything else (e.g. N).
/// With this coding the complement of code b is 3 - b.
#[inline]
pub fn base_code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Encode a sequence as a 2-bit packed integer. None if it has a non-ACGT base.
pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    let mut x = 0u64;
    for &b in seq {
        x = (x << 2) | base_code(b)?;
    }
    Some(x)
}

/// Reverse complement of a packed k-mer of length `k`.
#[inline]
pub fn revcomp_packed(mut x: u64, k: usize) -> u64 {
    let mut rc = 0u64;
    for _ in 0..k {
        rc = (rc << 2) | (3 - (x & 3)); // take the last base, complement it, append
        x >>= 2;
    }
    rc
}

/// Canonical form of a packed k-mer: the smaller of itself and its reverse complement.
pub fn canonical(x: u64, k: usize) -> u64 {
    x.min(revcomp_packed(x, k))
}

/// Canonical packed form of a DNA string (None if it contains non-ACGT bases).
pub fn canonical_of(seq: &[u8]) -> Option<u64> {
    encode_kmer(seq).map(|x| canonical(x, seq.len()))
}

/// The set of k-mers being searched for. Each distinct (k, canonical k-mer) gets a dense
/// index 0..len, so counts can live in a plain `Vec<u64>`.
///
/// Almost always there is a single k (2*flank+1). More than one only occurs when a
/// chromosome is shorter than the full window ("truncated" variants).
pub struct KmerIndex {
    /// k values in use, parallel to `maps`
    ks: Vec<usize>,
    maps: Vec<FastMap>,
    /// total number of distinct k-mers over all k
    n: usize,
}

impl KmerIndex {
    /// Empty set.
    pub fn new() -> Self {
        KmerIndex { ks: Vec::new(), maps: Vec::new(), n: 0 }
    }

    /// Add a canonical k-mer of length `k` (1..=31) and return its index
    /// (the existing index if already present).
    pub fn insert(&mut self, canon: u64, k: usize) -> usize {
        assert!((1..=31).contains(&k), "k must be between 1 and 31");
        let slot = match self.ks.iter().position(|&x| x == k) {
            Some(s) => s,
            None => {
                self.ks.push(k);
                self.maps.push(FastMap::default());
                self.ks.len() - 1
            }
        };
        let next = self.n;
        let i = *self.maps[slot].entry(canon).or_insert(next);
        if i == next {
            self.n += 1;
        }
        i
    }

    /// Index of a canonical k-mer of length `k`, if it is in the set.
    pub fn get(&self, canon: u64, k: usize) -> Option<usize> {
        let slot = self.ks.iter().position(|&x| x == k)?;
        self.maps[slot].get(&canon).copied()
    }

    /// Number of distinct k-mers.
    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }
}

impl Default for KmerIndex {
    fn default() -> Self {
        Self::new()
    }
}

/// Count every search k-mer occurring in one read, adding to `counts`.
///
/// A rolling window keeps the forward k-mer and its reverse complement up to date as
/// each base arrives, so each position costs a couple of shifts and one hash lookup.
/// With several k values the read is scanned once per k.
pub fn count_read(read: &[u8], index: &KmerIndex, counts: &mut [u64]) {
    for (&k, map) in index.ks.iter().zip(&index.maps) {
        let mask: u64 = (1u64 << (2 * k)) - 1;
        let rc_shift = 2 * (k - 1);
        let mut fwd = 0u64; // last k bases, forward strand
        let mut rc = 0u64; // reverse complement of the same k bases
        let mut valid = 0usize; // number of consecutive valid bases ending here

        for &b in read {
            match base_code(b) {
                Some(c) => {
                    fwd = ((fwd << 2) | c) & mask;
                    rc = (rc >> 2) | ((3 - c) << rc_shift);
                    valid += 1;
                    if valid >= k {
                        // canonical = the smaller strand representation
                        if let Some(&i) = map.get(&fwd.min(rc)) {
                            counts[i] += 1;
                        }
                    }
                }
                // N or other symbol: the window is broken, start again
                None => valid = 0,
            }
        }
    }
}

/// Open a FASTQ file, transparently handling gzip (detected by the magic bytes).
fn open_reads(path: &str) -> Result<Box<dyn BufRead + Send>, String> {
    let f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let mut br = BufReader::with_capacity(1 << 20, f);
    let is_gz = {
        let head = br.fill_buf().map_err(|e| format!("read error in {path}: {e}"))?;
        head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
    };
    if is_gz {
        // MultiGzDecoder also copes with files made of several concatenated gzip members
        Ok(Box::new(BufReader::with_capacity(1 << 20, MultiGzDecoder::new(br))))
    } else {
        Ok(Box::new(br))
    }
}

/// Count search k-mers across all reads of a FASTQ(.gz) file.
///
/// One thread decompresses and extracts the sequence lines (every 4th line of a FASTQ);
/// `threads` worker threads count k-mers in batches. Returns one count per index entry.
pub fn count_fastq(path: &str, index: &KmerIndex, threads: usize) -> Result<Vec<u64>, String> {
    let mut reader = open_reads(path)?;
    let n = index.len();
    let threads = threads.max(1);

    // batches of reads (sequences separated by '\n') are passed to the workers.
    // The channel is bounded so the reader cannot run far ahead of the counters.
    let (tx, rx) = sync_channel::<Vec<u8>>(threads * 2);
    let rx = Arc::new(Mutex::new(rx));

    thread::scope(|s| -> Result<Vec<u64>, String> {
        let workers: Vec<_> = (0..threads)
            .map(|_| {
                let rx = Arc::clone(&rx);
                s.spawn(move || {
                    let mut local = vec![0u64; n];
                    loop {
                        // take the lock only to receive, not while counting
                        let batch = match rx.lock().unwrap().recv() {
                            Ok(b) => b,
                            Err(_) => break, // reader finished
                        };
                        for read in batch.split(|&b| b == b'\n') {
                            count_read(read, index, &mut local);
                        }
                    }
                    local
                })
            })
            .collect();

        // reader (this thread): pick out the sequence line of each 4-line record
        let mut line = Vec::new();
        let mut batch: Vec<u8> = Vec::with_capacity(1 << 20);
        let mut line_no = 0u64;
        let mut read_err: Option<String> = None;
        loop {
            line.clear();
            match reader.read_until(b'\n', &mut line) {
                Ok(0) => break,
                Ok(_) => {
                    if line_no % 4 == 1 {
                        // strip the trailing newline / carriage return
                        while matches!(line.last(), Some(b'\n') | Some(b'\r')) {
                            line.pop();
                        }
                        batch.extend_from_slice(&line);
                        batch.push(b'\n');
                        if batch.len() >= 1 << 20 {
                            let full = std::mem::replace(&mut batch, Vec::with_capacity(1 << 20));
                            if tx.send(full).is_err() {
                                break;
                            }
                        }
                    }
                    line_no += 1;
                }
                Err(e) => {
                    read_err = Some(format!("read error in {path}: {e}"));
                    break;
                }
            }
        }
        if !batch.is_empty() {
            let _ = tx.send(batch);
        }
        // closing the sender tells the workers there is no more data
        drop(tx);

        // add up the per-thread counts
        let mut total = vec![0u64; n];
        for w in workers {
            for (t, c) in total.iter_mut().zip(w.join().unwrap()) {
                *t += c;
            }
        }
        match read_err {
            Some(e) => Err(e),
            None => Ok(total),
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::write::GzEncoder;
    use flate2::Compression;
    use std::io::Write;

    /// Index holding the canonical form of each given sequence.
    fn index_of(seqs: &[&str]) -> KmerIndex {
        let mut idx = KmerIndex::new();
        for s in seqs {
            idx.insert(canonical_of(s.as_bytes()).unwrap(), s.len());
        }
        idx
    }

    fn count(read: &str, idx: &KmerIndex) -> Vec<u64> {
        let mut c = vec![0u64; idx.len()];
        count_read(read.as_bytes(), idx, &mut c);
        c
    }

    #[test]
    fn revcomp_and_canonical() {
        // revcomp(ACGTT) = AACGT
        let x = encode_kmer(b"ACGTT").unwrap();
        assert_eq!(revcomp_packed(x, 5), encode_kmer(b"AACGT").unwrap());
        // a k-mer and its reverse complement share one canonical form
        assert_eq!(canonical_of(b"ACGTT"), canonical_of(b"AACGT"));
    }

    #[test]
    fn matches_either_strand_exactly() {
        let idx = index_of(&["ACGTT"]);
        assert_eq!(count("GGACGTTGG", &idx), vec![1]); // forward strand
        assert_eq!(count("GGAACGTGG", &idx), vec![1]); // reverse complement
        assert_eq!(count("GGACGATGG", &idx), vec![0]); // one mismatch: not counted
    }

    #[test]
    fn counts_repeats_and_read_ends() {
        let idx = index_of(&["ACGTT"]);
        assert_eq!(count("ACGTTACGTT", &idx), vec![2]); // twice, one at each end
        assert_eq!(count("ACGT", &idx), vec![0]); // read shorter than k
    }

    #[test]
    fn n_breaks_the_window() {
        let idx = index_of(&["ACGTT"]);
        assert_eq!(count("ACGNTT", &idx), vec![0]);
        assert_eq!(count("ACGTTNACGTT", &idx), vec![2]);
    }

    #[test]
    fn several_k_values_are_counted_separately() {
        let mut idx = KmerIndex::new();
        idx.insert(canonical_of(b"ACGTT").unwrap(), 5);
        idx.insert(canonical_of(b"ACGTTGC").unwrap(), 7);
        assert_eq!(idx.len(), 2);
        let mut c = vec![0u64; 2];
        count_read(b"GGACGTTGCGG", &idx, &mut c);
        assert_eq!(c, vec![1, 1]);
    }

    #[test]
    fn palindrome_counted_once_per_occurrence() {
        // ACGT is its own reverse complement
        let idx = index_of(&["ACGT"]);
        assert_eq!(count("TTACGTTT", &idx), vec![1]);
    }

    #[test]
    fn lower_case_reads_match() {
        let idx = index_of(&["ACGTT"]);
        assert_eq!(count("ggacgttgg", &idx), vec![1]);
    }

    /// Write a small FASTQ (plain and gzipped) and check counting works the same for both,
    /// and for any number of threads.
    #[test]
    fn fastq_plain_and_gz_any_threads() {
        let dir = std::env::temp_dir().join(format!("searchsnp_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let mut fq = String::new();
        for i in 0..5000 {
            // alternate forward and reverse-complement reads of the target
            let s = if i % 2 == 0 { "GGACGTTGG" } else { "GGAACGTGG" };
            fq.push_str(&format!("@r{i}\n{s}\n+\nIIIIIIIII\n"));
        }
        let plain = dir.join("r.fastq");
        std::fs::write(&plain, &fq).unwrap();
        let gz = dir.join("r.fastq.gz");
        let mut enc = GzEncoder::new(File::create(&gz).unwrap(), Compression::default());
        enc.write_all(fq.as_bytes()).unwrap();
        enc.finish().unwrap();

        let idx = index_of(&["ACGTT", "TTTTT"]);
        for threads in [1, 4] {
            for p in [&plain, &gz] {
                let c = count_fastq(p.to_str().unwrap(), &idx, threads).unwrap();
                assert_eq!(c, vec![5000, 0]);
            }
        }
        std::fs::remove_dir_all(&dir).ok();
    }
}
