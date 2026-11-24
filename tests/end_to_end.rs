//! End-to-end test: simulate a genome, variants (SNP, deletion, insertion, multi-allelic)
//! and reads with KNOWN numbers of ref/alt reads, run the real binary, and check the table.

use std::fs;
use std::path::PathBuf;
use std::process::Command;

/// Tiny deterministic pseudo-random generator so the test needs no extra crates.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        self.0 >> 33
    }
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            _ => 'C',
        })
        .collect()
}

/// FASTQ text with `n` copies of each read, every second one reverse-complemented
/// (reads come from either strand).
fn fastq(reads: &[(String, usize)]) -> String {
    let mut out = String::new();
    let mut id = 0;
    for (seq, n) in reads {
        for i in 0..*n {
            let s = if i % 2 == 0 { seq.clone() } else { revcomp(seq) };
            out.push_str(&format!("@r{id}\n{s}\n+\n{}\n", "I".repeat(s.len())));
            id += 1;
        }
    }
    out
}

fn tmpdir(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("searchsnp_e2e_{name}_{}", std::process::id()));
    fs::create_dir_all(&d).unwrap();
    d
}

/// 1 chromosome of 3000 random bases.
fn genome() -> String {
    let mut g = Lcg(42);
    (0..3000).map(|_| ['A', 'C', 'G', 'T'][(g.next() % 4) as usize]).collect()
}

/// A read: 20 reference bases before POS, the allele, then 20 bases after the REF allele.
fn read_with(genome: &str, pos: usize, ref_len: usize, allele: &str) -> String {
    format!("{}{}{}", &genome[pos - 21..pos - 1], allele, &genome[pos - 1 + ref_len..pos - 1 + ref_len + 20])
}

#[test]
fn known_read_counts_give_known_proportions() {
    let g = genome();
    let dir = tmpdir("counts");

    // variants: (pos, REF, ALT, number of ref reads, number of alt reads)
    let snp_pos = 500;
    let del_pos = 1000;
    let ins_pos = 1500;
    let multi_pos = 2000;
    let b = |p: usize, n: usize| g[p - 1..p - 1 + n].to_string();
    let snp_ref = b(snp_pos, 1);
    let snp_alt = if snp_ref == "A" { "C" } else { "A" }.to_string();
    let del_ref = b(del_pos, 4); // delete the 3 bases after the anchor
    let del_alt = b(del_pos, 1);
    let ins_ref = b(ins_pos, 1);
    let ins_alt = format!("{ins_ref}TTGA"); // insert 4 bases
    let multi_ref = b(multi_pos, 1);
    let alts: Vec<String> = ["A", "C", "G", "T"].iter().filter(|x| **x != multi_ref).map(|s| s.to_string()).collect();
    let (m1, m2) = (alts[0].clone(), alts[1].clone());

    // reference FASTA and VCF
    fs::write(dir.join("ref.fa"), format!(">chr1 test\n{g}\n")).unwrap();
    let vcf = format!(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\n\
         chr1\t{snp_pos}\t.\t{snp_ref}\t{snp_alt}\n\
         chr1\t{del_pos}\t.\t{del_ref}\t{del_alt}\n\
         chr1\t{ins_pos}\t.\t{ins_ref}\t{ins_alt}\n\
         chr1\t{multi_pos}\t.\t{multi_ref}\t{m1},{m2}\n"
    );
    fs::write(dir.join("v.vcf"), vcf).unwrap();

    // reads with known counts
    let reads = vec![
        (read_with(&g, snp_pos, 1, &snp_ref), 30),
        (read_with(&g, snp_pos, 1, &snp_alt), 10),
        (read_with(&g, del_pos, 4, &del_ref), 20),
        (read_with(&g, del_pos, 4, &del_alt), 20),
        (read_with(&g, ins_pos, 1, &ins_ref), 8),
        (read_with(&g, ins_pos, 1, &ins_alt), 32),
        (read_with(&g, multi_pos, 1, &multi_ref), 10),
        (read_with(&g, multi_pos, 1, &m1), 10), // m2 has no reads: a "false positive"
    ];
    fs::write(dir.join("fwd.fastq"), fastq(&reads)).unwrap();
    fs::write(dir.join("rev.fastq"), fastq(&reads)).unwrap();

    let out = Command::new(env!("CARGO_BIN_EXE_searchSNP"))
        .args(["--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "fwd.fastq", "--rev", "rev.fastq"])
        .args(["--flank", "11", "--threads", "2"])
        .current_dir(&dir)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));

    let table = fs::read_to_string(dir.join("searchSNP_summary.csv")).unwrap();
    let rows: Vec<Vec<&str>> = table.lines().skip(1).map(|l| l.split(',').collect()).collect();
    assert_eq!(rows.len(), 5, "4 VCF rows, one of them multi-allelic -> 5 rows");

    // columns: 0 chrom,1 pos,2 ref,3 alt,4 ref_matches,5 multiallelic,6 n_alts,7 window_bp,
    
    let check = |row: &Vec<&str>, r: &str, a: &str, prop: &str| {
        assert_eq!((row[4], row[5]), (r, a), "fwd counts of {row:?}");
        assert_eq!((row[8], row[9]), (r, a), "rev counts of {row:?}");
        assert_eq!(row[7], prop);
        assert_eq!((row[12], row[17]), ("23", "full"));
        assert_eq!(row[13], "true");
    };
    check(&rows[0], "30", "10", "0.250000"); // SNP
    check(&rows[1], "20", "20", "0.500000"); // deletion
    check(&rows[2], "8", "32", "0.800000"); // insertion
    check(&rows[3], "10", "10", "0.500000"); // multi-allelic, first alt
    check(&rows[4], "10", "0", "0.000000"); // multi-allelic, second alt: no support
    assert_eq!((rows[3][14], rows[3][15]), ("true", "2"));
    assert_eq!((rows[0][14], rows[0][15]), ("false", "1"));
    // output stays in VCF order
    assert_eq!(rows.iter().map(|r| r[1]).collect::<Vec<_>>(), ["500", "1000", "1500", "2000", "2000"]);
    fs::remove_dir_all(&dir).ok();
}

#[test]
fn wrong_reference_is_flagged_and_strict_stops() {
    let g = genome();
    let dir = tmpdir("strict");
    fs::write(dir.join("ref.fa"), format!(">chr1\n{g}\n")).unwrap();
    // REF deliberately wrong (reference base + a different base)
    let real = &g[599..600];
    let wrong = if real == "A" { "C" } else { "A" };
    fs::write(
        dir.join("v.vcf"),
        format!("#CHROM\tPOS\tID\tREF\tALT\nchr1\t600\t.\t{wrong}\tT\n"),
    )
    .unwrap();
    let r = read_with(&g, 600, 1, "T");
    fs::write(dir.join("r.fastq"), fastq(&[(r, 4)])).unwrap();

    let run = |extra: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_searchSNP"))
            .args(["--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "r.fastq", "--rev", "r.fastq"])
            .args(extra)
            .current_dir(&dir)
            .output()
            .unwrap()
    };
    // default: warns, carries on, and ref_matches is false
    let ok = run(&[]);
    assert!(ok.status.success());
    assert!(String::from_utf8_lossy(&ok.stderr).contains("REF mismatch"));
    let table = fs::read_to_string(dir.join("searchSNP_summary.csv")).unwrap();
    assert!(table.lines().nth(1).unwrap().contains(",false,"));
    // --strict: fails
    assert!(!run(&["--strict"]).status.success());
    // flank outside 1..=15 is rejected
    assert!(!run(&["--flank", "16"]).status.success());
    fs::remove_dir_all(&dir).ok();
}

/// A variant 3 bases from the contig start is searched with a shifted, full-length window,
/// and a contig shorter than the minimum window is reported as skipped (not as "0 reads").
#[test]
fn edge_variants_are_shifted_or_skipped_and_reported() {
    let g = genome();
    let dir = tmpdir("edge");
    fs::write(dir.join("ref.fa"), format!(">chr1\n{g}\n>tiny\nACGTACGTAC\n")).unwrap();
    let snp_ref = &g[2..3];
    let snp_alt = if snp_ref == "A" { "C" } else { "A" };
    fs::write(
        dir.join("v.vcf"),
        format!("#CHROM\tPOS\tID\tREF\tALT\nchr1\t3\t.\t{snp_ref}\t{snp_alt}\ntiny\t5\t.\tA\tG\n"),
    )
    .unwrap();
    // reads: 12 reference bases from base 1 (the shifted window is bases 1..23), plus alt
    let rd = |allele: &str| format!("{}{}{}", &g[0..2], allele, &g[3..40]);
    fs::write(dir.join("r.fastq"), fastq(&[(rd(snp_ref), 6), (rd(snp_alt), 2)])).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_searchSNP"))
        .args(["--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "r.fastq", "--rev", "r.fastq"])
        .current_dir(&dir)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let table = fs::read_to_string(dir.join("searchSNP_summary.csv")).unwrap();
    let rows: Vec<Vec<&str>> = table.lines().skip(1).map(|l| l.split(',').collect()).collect();
    assert_eq!((rows[0][12], rows[0][17]), ("23", "shifted"));
    assert_eq!((rows[0][4], rows[0][5]), ("6", "2"));
    assert_eq!((rows[1][12], rows[1][17], rows[1][16]), ("0", "skipped_too_short", "0"));
    assert_eq!((rows[1][4], rows[1][5], rows[1][7]), ("0", "0", "0.000000"));
    fs::remove_dir_all(&dir).ok();
}

/// sequence_counts.csv lists every searched sequence with its counts, and a placeholder
/// row for a variant skipped by a rule.
#[test]
fn sequence_counts_table_has_rows_and_placeholders() {
    let g = genome();
    let dir = tmpdir("seqcounts");
    fs::write(dir.join("ref.fa"), format!(">chr1\n{g}\n>tiny\nACGTACGTAC\n")).unwrap();
    let snp_ref = &g[999..1000];
    let snp_alt = if snp_ref == "A" { "C" } else { "A" };
    fs::write(
        dir.join("v.vcf"),
        format!("#CHROM\tPOS\tID\tREF\tALT\nchr1\t1000\t.\t{snp_ref}\t{snp_alt}\ntiny\t5\t.\tA\tG\n"),
    )
    .unwrap();
    let rd = |a: &str| read_with(&g, 1000, 1, a);
    fs::write(dir.join("fwd.fastq"), fastq(&[(rd(snp_ref), 7), (rd(snp_alt), 3)])).unwrap();
    fs::write(dir.join("rev.fastq"), fastq(&[(rd(snp_ref), 5)])).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_searchSNP"))
        .args(["--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "fwd.fastq", "--rev", "rev.fastq"])
        .current_dir(&dir)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let t = fs::read_to_string(dir.join("searchSNP_sequence_counts.csv")).unwrap();
    let rows: Vec<Vec<&str>> = t.lines().skip(1).map(|l| l.split(',').collect()).collect();
    assert_eq!(rows.len(), 3); // ref, alt, and the skipped placeholder
    assert_eq!((rows[0][4], rows[0][6], rows[0][7]), ("ref", "7", "5"));
    assert_eq!((rows[1][4], rows[1][6], rows[1][7]), ("alt", "3", "0"));
    assert_eq!((rows[2][5], rows[2][9], rows[2][6], rows[2][7]), ("", "skipped_too_short", "", ""));
    assert!(!dir.join("searchSNP_fwd_sequences.fasta").exists());
    fs::remove_dir_all(&dir).ok();
}

/// A long insertion next to a focal deletion fills the window, so the insertion haplotype
/// gives the same k-mer under both labels. It must be flagged per sequence and left out of
/// the deletion's counts, while clean reads still count.
#[test]
fn ambiguous_sequence_is_flagged_and_excluded() {
    let g = genome();
    let dir = tmpdir("ambig");
    fs::write(dir.join("ref.fa"), format!(">chr1\n{g}\n")).unwrap();
    let b = |p: usize, n: usize| g[p - 1..p - 1 + n].to_string();
    let ins_ref = b(1000, 1);
    let ins_alt = format!("{ins_ref}{}", "ACGTTGCAAGGTCCATTGACCTAGGA"); // longer than the window
    let del_ref = b(1003, 3);
    let del_alt = b(1003, 1);
    fs::write(
        dir.join("v.vcf"),
        format!("#CHROM\tPOS\tID\tREF\tALT\nchr1\t1003\t.\t{del_ref}\t{del_alt}\nchr1\t1000\t.\t{ins_ref}\t{ins_alt}\n"),
    )
    .unwrap();
    // reads: 6 clean ref, 4 clean deletion, 50 carrying the insertion (cannot tell ref/alt)
    let ins_read = format!("{}{}{}", &g[979..999], ins_alt, &g[1000..1020]);
    let reads = vec![
        (read_with(&g, 1003, 3, &del_ref), 6),
        (read_with(&g, 1003, 3, &del_alt), 4),
        (ins_read, 50),
    ];
    fs::write(dir.join("r.fastq"), fastq(&reads)).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_searchSNP"))
        .args(["--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "r.fastq", "--rev", "r.fastq"])
        .current_dir(&dir)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    // summary: the deletion's counts exclude the 50 insertion reads; 1 ambiguous sequence
    let t = fs::read_to_string(dir.join("searchSNP_summary.csv")).unwrap();
    let row: Vec<&str> = t.lines().nth(1).unwrap().split(',').collect();
    assert_eq!((row[4], row[5], row[19]), ("6", "4", "1"));
    // per-sequence table: the ambiguous sequence is flagged and shows its real count
    let t = fs::read_to_string(dir.join("searchSNP_sequence_counts.csv")).unwrap();
    let rows: Vec<Vec<&str>> = t.lines().skip(1).map(|l| l.split(',').collect()).collect();
    let del_rows: Vec<_> = rows.iter().filter(|r| r[1] == "1003").collect();
    // the same string appears once as ref and once as alt: both rows are flagged
    assert_eq!(del_rows.iter().filter(|r| r[11] == "true").count(), 2);
    assert!(del_rows.iter().any(|r| r[11] == "true" && r[6] == "50"));
    fs::remove_dir_all(&dir).ok();
}

/// summary `shared_in_reference`: which of the ref/alt search sequences also occur
/// somewhere else in the reference genome (`ref`, `alt`, `both` or `none`).
#[test]
fn shared_in_reference_labels_ref_alt_both_none() {
    let g0 = genome().into_bytes();
    let other_base = |b: u8| if b == b'A' { b'C' } else { b'A' };
    // 23 bp window (flank 11) around a 1-based position, focal base replaced by `allele`
    let window = |g: &[u8], pos: usize, allele: u8| {
        let mut w = g[pos - 12..pos + 11].to_vec();
        w[11] = allele;
        w
    };

    let (p_alt, p_both, p_none) = (1000usize, 2000usize, 500usize);
    let (r_alt, r_both, r_none) = (g0[p_alt - 1], g0[p_both - 1], g0[p_none - 1]);
    let (a_alt, a_both, a_none) = (other_base(r_alt), other_base(r_both), other_base(r_none));

    let w_alt = window(&g0, p_alt, a_alt);
    let w_both_ref = window(&g0, p_both, r_both);
    let w_both_alt = window(&g0, p_both, a_both);

    // variant 1: only its ALT window occurs elsewhere; variant 2: both windows do
    let mut g = g0;
    g[2499..2522].copy_from_slice(&w_alt);
    g[2749..2772].copy_from_slice(&w_both_ref);
    g[2899..2922].copy_from_slice(&w_both_alt);

    let dir = tmpdir("shared_ref");
    fs::write(dir.join("ref.fa"), format!(">chr1\n{}\n", String::from_utf8(g).unwrap())).unwrap();
    fs::write(
        dir.join("v.vcf"),
        format!(
            "#CHROM\tPOS\tID\tREF\tALT\n\
             chr1\t{p_alt}\t.\t{}\t{}\n\
             chr1\t{p_both}\t.\t{}\t{}\n\
             chr1\t{p_none}\t.\t{}\t{}\n",
            r_alt as char,
            a_alt as char,
            r_both as char,
            a_both as char,
            r_none as char,
            a_none as char
        ),
    )
    .unwrap();
    // the counts are irrelevant here, so the read files can be empty
    fs::write(dir.join("fwd.fastq"), "").unwrap();
    fs::write(dir.join("rev.fastq"), "").unwrap();

    // run and return the shared_in_reference column (the last one), in VCF order
    let column = |extra: &[&str], prefix: &str| -> Vec<String> {
        let out = Command::new(env!("CARGO_BIN_EXE_searchSNP"))
            .args([
                "--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "fwd.fastq", "--rev",
                "rev.fastq", "--prefix", prefix,
            ])
            .args(extra)
            .current_dir(&dir)
            .output()
            .unwrap();
        assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
        let t = fs::read_to_string(dir.join(format!("{prefix}_summary.csv"))).unwrap();
        t.lines().skip(1).map(|l| l.split(',').last().unwrap().to_string()).collect()
    };

    assert_eq!(column(&[], "checked"), ["alt", "both", "none"]);
    // nothing is checked with --no-ref-check, so the column is empty
    assert_eq!(column(&["--no-ref-check"], "unchecked"), ["", "", ""]);
    fs::remove_dir_all(&dir).ok();
}

/// summary: a variant skipped by a rule has no sequences, so `shared_in_reference` is empty
/// (not `none`, which would claim the reference had been checked).
#[test]
fn skipped_variant_has_empty_shared_in_reference() {
    let g = genome();
    let dir = tmpdir("skipped_shared");
    fs::write(dir.join("ref.fa"), format!(">chr1\n{g}\n>tiny\nACGTACGTAC\n")).unwrap();
    fs::write(dir.join("v.vcf"), "#CHROM\tPOS\tID\tREF\tALT\nchr1\t1500\t.\tA\tC\ntiny\t5\t.\tA\tG\n").unwrap();
    fs::write(dir.join("fwd.fastq"), "").unwrap();
    fs::write(dir.join("rev.fastq"), "").unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_searchSNP"))
        .args(["--vcf", "v.vcf", "--reference", "ref.fa", "--fwd", "fwd.fastq", "--rev", "rev.fastq"])
        .current_dir(&dir)
        .output()
        .unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    let t = fs::read_to_string(dir.join("searchSNP_summary.csv")).unwrap();
    let rows: Vec<Vec<&str>> = t.lines().skip(1).map(|l| l.split(',').collect()).collect();
    assert_eq!(rows[1][17], "skipped_too_short");
    assert_eq!(rows[1][20], "");
    fs::remove_dir_all(&dir).ok();
}
