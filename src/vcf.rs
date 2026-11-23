//! Reading variants from a VCF file.
//!
//! Each output `Variant` is one (position, ref, alt) triple. A VCF row with several
//! comma-separated ALT alleles is split into one `Variant` per allele, in order.

use crate::reference::Reference;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// A single variant (one ALT allele of one VCF row).
#[derive(Clone, Debug)]
pub struct Variant {
    /// Index of the chromosome in the `Reference`.
    pub chrom: usize,
    /// 1-based position of the first base of REF (the anchor base for indels).
    pub pos: usize,
    /// REF allele (upper case).
    pub r: Vec<u8>,
    /// ALT allele (upper case).
    pub a: Vec<u8>,
    /// Does the VCF REF equal the reference FASTA at POS? (case-insensitive)
    pub ref_matches: bool,
}

/// Result of reading a VCF.
pub struct VcfData {
    pub variants: Vec<Variant>,
    /// Number of variants whose REF did not match the reference.
    pub n_ref_mismatch: usize,
}

/// Parse a VCF. Fails if a chromosome is not in the reference, or the file is malformed.
/// REF mismatches are not errors here; they are flagged and counted (the caller decides).
pub fn read_vcf(path: &str, reference: &Reference) -> Result<VcfData, String> {
    let f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
    let mut variants = Vec::new();
    let mut n_ref_mismatch = 0;

    for line in BufReader::new(f).lines() {
        let line = line.map_err(|e| format!("read error in {path}: {e}"))?;
        // skip header / meta lines and blanks
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            return Err(format!("malformed VCF line (need >= 5 columns): {line}"));
        }
        // the chromosome must exist in the reference: wrong FASTA is the commonest mix-up
        let chrom = reference.index_of(f[0]).ok_or_else(|| {
            format!(
                "chromosome '{}' in the VCF is not in the reference FASTA (wrong reference file?)",
                f[0]
            )
        })?;
        let pos: usize = f[1].parse().map_err(|_| format!("bad POS '{}' in VCF", f[1]))?;
        if pos == 0 {
            return Err("POS must be >= 1".to_string());
        }
        let r = f[3].to_ascii_uppercase().into_bytes();

        // check REF against the FASTA once per row (same REF for all its ALTs)
        let seq = &reference.chroms[chrom].seq;
        let ref_matches = pos - 1 + r.len() <= seq.len() && seq[pos - 1..pos - 1 + r.len()] == r[..];

        let alts: Vec<&str> = f[4].split(',').collect();
        for alt in &alts {
            // symbolic / missing alleles cannot be turned into a sequence
            if alt.starts_with('<') || alt.contains('*') || *alt == "." {
                eprintln!("warning: skipping unsupported ALT '{alt}' at {}:{pos}", f[0]);
                continue;
            }
            if !ref_matches {
                n_ref_mismatch += 1;
                eprintln!(
                    "warning: REF mismatch at {}:{pos}: VCF has {}, reference has {}",
                    f[0],
                    String::from_utf8_lossy(&r),
                    String::from_utf8_lossy(
                        &seq[(pos - 1).min(seq.len())..(pos - 1 + r.len()).min(seq.len())]
                    )
                );
            }
            variants.push(Variant {
                chrom,
                pos,
                r: r.clone(),
                a: alt.to_ascii_uppercase().into_bytes(),
                ref_matches,
            });
        }
    }
    Ok(VcfData { variants, n_ref_mismatch })
}
