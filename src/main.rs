//! searchSNP command line program.
//!
//! Given a VCF of called variants, the reference genome it was called against, and the raw
//! reads (forward and reverse files), count exact matches of the ref and alt sequences around
//! every variant in the reads and write a per-variant summary table.

use clap::Parser;
use search_snp::filter::filter_pairs;
use search_snp::kmers::{count_fastq, KmerIndex};
use search_snp::reference::Reference;
use search_snp::refcheck::annotate;
use search_snp::sequences::build_all;
use search_snp::summary::{index_sequences, write_sequence_counts, write_summary};
use search_snp::vcf::read_vcf;
use search_snp::MAX_FLANK;
use std::fs::{self, File};
use std::io::BufWriter;
use std::process;
use std::thread;
use std::time::Instant;

#[derive(Parser)]
#[command(
    version,
    about = "Check whether called SNPs/indels are supported by exact k-mer matches in raw reads",
    long_about = "For every row of a VCF, builds the reference and alternate sequences around the \
variant (flank bases before the position, then the allele, then reference, always 2*flank+1 \
long) and counts exact matches of them (either strand) in the forward and reverse reads. Writes one \
summary row per VCF row."
)]
struct Args {
    /// VCF file of variants (e.g. from breseq)
    #[arg(long)]
    vcf: String,

    /// Reference FASTA the VCF was called against
    #[arg(long)]
    reference: String,

    /// Forward reads (e.g. the R1 file), FASTQ or FASTQ.gz
    #[arg(long)]
    fwd: String,

    /// Reverse reads (e.g. the R2 file), FASTQ or FASTQ.gz
    #[arg(long)]
    rev: String,

    /// Bases before (and after) the variant; sequences are 2*flank+1 long (max 15)
    #[arg(long, default_value_t = 11, value_parser = clap::value_parser!(u8).range(1..=MAX_FLANK as i64))]
    flank: u8,

    /// K-mer counting threads, split across the two read files (default: all
    /// available cores); one reader thread per file decompresses it
    #[arg(long)]
    threads: Option<usize>,

    /// Directory for the output files (created if missing)
    #[arg(long, default_value = ".")]
    outdir: String,

    /// Prefix for the output file names
    #[arg(long, default_value = "searchSNP")]
    prefix: String,

    /// Stop with an error if any VCF REF does not match the reference FASTA
    #[arg(long)]
    strict: bool,

    /// Skip the check of whether the search sequences also occur elsewhere in the reference
    #[arg(long)]
    no_ref_check: bool,

    /// Forward reads of the ancestor sample (optional, recommended). Haplotype pairs whose
    /// ALT sequence is found in these reads are left out of <prefix>_summary_filtered.csv
    #[arg(long, requires = "ancestor_rev")]
    ancestor_fwd: Option<String>,

    /// Reverse reads of the ancestor sample (required with --ancestor-fwd)
    #[arg(long, requires = "ancestor_fwd")]
    ancestor_rev: Option<String>,

    /// Ancestor ALT reads (forward + reverse) needed to drop a ref/alt haplotype pair
    #[arg(long, default_value_t = 1, requires = "ancestor_fwd")]
    min_alt_count: u64,

    /// Skip variants that would need more haplotype combinations than this
    #[arg(long, default_value_t = 1024)]
    max_combos: u128,
}

/// Print an error and exit with a non-zero status.
fn die(msg: &str) -> ! {
    eprintln!("error: {msg}");
    process::exit(1);
}

/// Count the k-mers of one forward/reverse read pair, splitting `threads` between the two
/// files. With a single worker there is nothing to share, so the files are counted one
/// after the other.
fn count_pair(fwd: &str, rev: &str, index: &KmerIndex, threads: usize) -> (Vec<u64>, Vec<u64>) {
    let fwd_threads = threads.div_ceil(2);
    let rev_threads = threads / 2;
    let (fwd_res, rev_res) = if rev_threads == 0 {
        (count_fastq(fwd, index, 1), count_fastq(rev, index, 1))
    } else {
        thread::scope(|s| {
            let h1 = s.spawn(|| count_fastq(fwd, index, fwd_threads));
            let h2 = s.spawn(|| count_fastq(rev, index, rev_threads));
            (h1.join().unwrap(), h2.join().unwrap())
        })
    };
    (fwd_res.unwrap_or_else(|e| die(&e)), rev_res.unwrap_or_else(|e| die(&e)))
}

fn main() {
    let args = Args::parse();
    let flank = args.flank as usize;
    let k = 2 * flank + 1;
    let threads = args
        .threads
        .unwrap_or_else(|| thread::available_parallelism().map(|n| n.get()).unwrap_or(1))
        .max(1);
    let t0 = Instant::now();

    // 1. reference and variants (with the REF-vs-FASTA check)
    let reference = Reference::from_fasta(&args.reference).unwrap_or_else(|e| die(&e));
    let vcf = read_vcf(&args.vcf, &reference).unwrap_or_else(|e| die(&e));
    let vars = &vcf.variants;
    eprintln!("{} variants, {} chromosomes", vars.len(), reference.chroms.len());
    if vcf.n_ref_mismatch > 0 {
        let msg = format!(
            "{} of {} variants have a REF that does not match the reference FASTA \
             (is this the reference the VCF was called against?)",
            vcf.n_ref_mismatch,
            vars.len()
        );
        if args.strict {
            die(&msg);
        }
        eprintln!("warning: {msg}");
    }

    // 2. search sequences for every variant
    let mut seqs = build_all(vars, &reference, flank, args.max_combos, threads);
    let (index, vk) = index_sequences(&seqs, vars);
    eprintln!("{} distinct k-mers (k={k}) to search for", index.len());

    // 2b. reference-occurrence check: does each search sequence also occur elsewhere?
    if !args.no_ref_check {
        let t_ref = Instant::now();
        let n_else = annotate(&mut seqs, vars, &reference);
        let n_searchable: usize = seqs.iter().map(|s| s.seqs.len()).sum();
        eprintln!(
            "reference occurrence check: {n_else} of {n_searchable} search sequences occur \
             elsewhere in the reference ({:.1}s)",
            t_ref.elapsed().as_secs_f64()
        );
    }

    let n_skipped = seqs.iter().filter(|s| s.seqs.is_empty()).count();
    let n_edge = seqs.iter().filter(|s| s.status != search_snp::sequences::SearchStatus::Full && !s.seqs.is_empty()).count();
    if n_edge + n_skipped > 0 {
        eprintln!("{n_edge} variants near a chromosome end (shifted/truncated window), {n_skipped} skipped (see search_status)");
    }

    // 3. count them in the forward and reverse reads (the two files are read at the same time).
    // `threads` is the number of worker threads doing the k-mer lookups; each file also has a
    // reader thread of its own that decompresses and parses it. With a single worker there is
    // nothing to share, so the two files are counted one after the other.
    let (counts_fwd, counts_rev) = count_pair(&args.fwd, &args.rev, &index, threads);

    // 3b. optional: count the ancestral reads with the same index, so a pair means the same
    // thing in both samples
    let ancestor = match (&args.ancestor_fwd, &args.ancestor_rev) {
        (Some(af), Some(ar)) => {
            let t_anc = Instant::now();
            let counts = count_pair(af, ar, &index, threads);
            eprintln!("ancestor reads counted ({:.1}s)", t_anc.elapsed().as_secs_f64());
            Some(counts)
        }
        _ => None,
    };

    // 3c. remove the haplotype pairs whose ALT is present in the ancestor, and rebuild the
    // per-variant k-mer lists from what remains
    let filtered = ancestor.as_ref().map(|(anc_fwd, anc_rev)| {
        let (fseqs, fkmers, status, stats) =
            filter_pairs(&seqs, vars, &index, anc_fwd, anc_rev, args.min_alt_count);
        eprintln!(
            "ancestor filter: {} of {} pairs removed; {} variants partial, {} with no pairs left",
            stats.dropped, stats.pairs, stats.partial, stats.all_dropped
        );
        (fseqs, fkmers, status)
    });

    // 4. write the outputs
    fs::create_dir_all(&args.outdir).unwrap_or_else(|e| die(&format!("cannot create {}: {e}", args.outdir)));
    let open = |suffix: &str| {
        let p = format!("{}/{}_{}", args.outdir, args.prefix, suffix);
        BufWriter::new(File::create(&p).unwrap_or_else(|e| die(&format!("cannot write {p}: {e}"))))
    };
    write_summary(open("summary.csv"), vars, &vk, &seqs, &reference, &counts_fwd, &counts_rev, !args.no_ref_check, None)
        .unwrap_or_else(|e| die(&format!("write error: {e}")));
    write_sequence_counts(
        open("sequence_counts.csv"),
        vars,
        &seqs,
        &vk,
        &index,
        &reference,
        &counts_fwd,
        &counts_rev,
        !args.no_ref_check,
        ancestor.as_ref().map(|(a, b)| (a.as_slice(), b.as_slice())),
    )
    .unwrap_or_else(|e| die(&format!("write error: {e}")));
    if let Some((fseqs, fkmers, status)) = &filtered {
        write_summary(
            open("summary_filtered.csv"),
            vars,
            fkmers,
            fseqs,
            &reference,
            &counts_fwd,
            &counts_rev,
            !args.no_ref_check,
            Some(status),
        )
        .unwrap_or_else(|e| die(&format!("write error: {e}")));
    }

    eprintln!("done in {:.1}s, wrote {}/{}_summary.csv", t0.elapsed().as_secs_f64(), args.outdir, args.prefix);
}
