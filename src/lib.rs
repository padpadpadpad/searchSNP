//! searchSNP library: the building blocks of the pipeline.
//!
//! Pipeline: read reference + VCF -> build search sequences around every variant
//! -> count exact (canonical) k-mer matches in forward/reverse reads -> summarise per VCF row.

pub mod filter;
pub mod kmers;
pub mod reference;
pub mod refcheck;
pub mod sequences;
pub mod summary;
pub mod vcf;

/// Largest allowed flank. k = 2*flank + 1 = 31 is the largest odd k that fits
/// in a u64 (2 bits per base => 32 bases max).
pub const MAX_FLANK: usize = 15;
