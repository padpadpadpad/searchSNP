//! Reading the reference genome from a FASTA file.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};

/// One chromosome / contig of the reference.
pub struct Chrom {
    /// Name = first word of the FASTA header (this is what the VCF CHROM column must match).
    pub name: String,
    /// Upper-case sequence bytes.
    pub seq: Vec<u8>,
}

/// The whole reference plus a name -> index lookup.
pub struct Reference {
    pub chroms: Vec<Chrom>,
    index: HashMap<String, usize>,
}

impl Reference {
    /// Build a reference from already-parsed chromosomes (used by tests).
    pub fn from_chroms(chroms: Vec<Chrom>) -> Self {
        let index = chroms
            .iter()
            .enumerate()
            .map(|(i, c)| (c.name.clone(), i))
            .collect();
        Reference { chroms, index }
    }

    /// Read a FASTA file. Sequences are upper-cased so soft-masked (lower-case)
    /// bases do not matter.
    pub fn from_fasta(path: &str) -> Result<Self, String> {
        let f = File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?;
        let mut chroms: Vec<Chrom> = Vec::new();
        for line in BufReader::new(f).lines() {
            let line = line.map_err(|e| format!("read error in {path}: {e}"))?;
            if let Some(h) = line.strip_prefix('>') {
                // start a new chromosome; the name stops at the first whitespace
                let name = h.split_whitespace().next().unwrap_or("").to_string();
                chroms.push(Chrom { name, seq: Vec::new() });
            } else if let Some(c) = chroms.last_mut() {
                // sequence line: append (upper-cased) to the current chromosome
                c.seq.extend(line.trim_end().bytes().map(|b| b.to_ascii_uppercase()));
            }
        }
        if chroms.is_empty() {
            return Err(format!("no sequences found in {path}"));
        }
        Ok(Self::from_chroms(chroms))
    }

    /// Index of a chromosome by name, if present.
    pub fn index_of(&self, name: &str) -> Option<usize> {
        self.index.get(name).copied()
    }
}
