# searchSNP

**searchSNP** checks whether called SNPs and indels in a vcf are real by looking for them directly in the raw reads that were used to generate the vcf.

**This is a development package and if you use it, we recommend checking the findings make sense for your usage (e.g. compare the proportions you get from searchSNP against those from your variant caller). We welcome any [Issues](https://github.com/padpadpadpad/searchSNP/issues) to improve the workflow and documentation.**

Calling SNPs in mixed populations is renowned for high rates of false positives. These can be hard to distinguish and filter. This tool attempts to overcome this by searching for short kmers of each reference and alternate sequence in the raw reads. For every row of a VCF from a variant caller (e.g. from [breseq](https://github.com/barricklab/breseq)) it builds the reference and alternate DNA sequence around the variant and counts **exact** matches of those sequences in the forward and reverse reads. From this, you can compare the presence/absence and proportion to that of your variant caller.

There are many filters that can be applied to attempt in an attempt to balance false positives and false negatives. Many of these are summarised and discussed in [Bush et al. 2021](https://doi.org/10.1099/mgen.0.000615). 

For population-level sequencing of experimental evolution studies, we prefer a stringent approach where we might lose some true positives (have some false negatives), to minimise false positives.

Consequently, we believe that **searchSNP** can give the following insight:

- If the predicted variant kmer **IS NOT** in any of the raw reads, then it is likely a false positive as it lacks direct read-level support. **We recommend removing predicted variants with no support in the raw reads.**
- If the predicted variant kmer **IS** in the reads of the ancestral/reference genome, the variant should be interpreted with caution. Detection of the predicted variant kmer in the ancestor indicates that the sequence was already present before the experiment, either due to pre-existing variation, homologous sequence, or discrepancies between the ancestral genome and the reference assembly. **For analyses focused on identifying novel mutations, we recommend excluding ALT/REF k-mer pairs where the ALT k-mer is detected in ancestral reads.**
- If the predicted variant kmer is only present on one strand (forward or reverse reads) then it could be a false positive. We recommend filtering out variants that are only present in one strand, as strand bias, where the alternate allele is disproportionately found on either the forward or reverse reads, increases the likelihood of a false-positive call [Guo et al 2012](https://doi.org/10.1186/1471-2164-13-666).

## Install

Prebuilt self-contained binaries for Windows, macOS and Linux are attached to every
[release](https://github.com/padpadpadpad/searchSNP/releases). They need no Rust or other software installed.

Download, unpack, and put `searchSNP` (`searchSNP.exe` on Windows) anywhere on your `PATH`. Or use a one-line installer:

**macOS and Linux** - installs to `~/.local/bin` (set `INSTALL_DIR` to change):

```sh
curl -fsSL https://raw.githubusercontent.com/padpadpadpad/searchSNP/main/install.sh | bash
```

**Windows (PowerShell)** - installs to `%LOCALAPPDATA%\searchSNP\bin` and adds it to your user `PATH`
(set `$env:INSTALL_DIR` to change); open a new terminal afterwards:

```powershell
irm https://raw.githubusercontent.com/padpadpadpad/searchSNP/main/install.ps1 | iex
```

**With Rust** (any platform):

```sh
cargo install --git https://github.com/padpadpadpad/searchSNP
```

**From source:**

```sh
git clone https://github.com/padpadpadpad/searchSNP
cd searchSNP
cargo build --release
# binary: target/release/searchSNP
```

**macOS note:** if you downloaded a binary with a browser rather than the installer, macOS may refuse to run it ("cannot be opened because the developer cannot be verified"). Right-click the file and choose Open, or run `xattr -d com.apple.quarantine searchSNP`.

**If `searchSNP` is not found:**

- In the folder you unpacked it into, run `./searchSNP` - Linux and macOS do not look in the current directory, so the name on its own is not enough.
- Straight after the installers, the binary is at `~/.local/bin/searchSNP`. Run that full path, or add the directory to your `PATH` and open a new terminal. Use `~/.zshrc` for zsh, `~/.bashrc` for bash:

```sh
# bash
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc

# zsh
echo 'export PATH="$HOME/.local/bin:$PATH"' >> ~/.zshrc && source ~/.zshrc
```

- Check the install with `searchSNP --version`.

## Licensing

This program is licensed under an MIT license.

## AI use

This program is a re-write of janky R code written by Daniel Padfield to create the search sequences and [jellyfish](https://github.com/gmarcais/jellyfish). Open Code with GitHub Copilot (Claude Sonnet 5.5) was used to rewrite all the steps in Rust. The free model GLM-5.3-Flash Free was used to help write validations.

## Idea

The idea behind **searchSNP** is to provide a way to validate the results of variant calling by directly searching for the variants in the raw sequencing data.

Let's say you have the a VCF row that looks like:

| #CHROM | POS | ID | REF | ALT | QUAL | FILTER | INFO |
|---------|-----|----|-----|-----|------|--------|-------------|
| contig1 | 111100 | . | T | C | 500 | PASS | AF=1.0000 |

This means we have a T -> C substitution at position 111100 on contig1, and it is at fixation (AF [alternate frequency] = 1).

If this is a real variant, then we would expect a short kmer of this to be present in the raw reads (typically between 150-250 base pairs for short reads). And the reference kmer should not be there. This is the simple idea of searchSNP.

![](diagram.png)

## General rules

- If the SNP has other SNPs up or downstream within the flanking region, **searchSNP** creates all possible sequences of the reference and alternate sequences. This scales with 2^n with n being the number of other variants in the region.
- If there are multiple allelic alternates (alt = A,T), then **searchSNP** will return two rows for the position, not one.
- Deletions stitch together either side of the deletion.
- Insertions take from the left (upstream flank) and the focal base pair is the first base pair before the insertion, followed by the start of the insertion which becomes the downstream flank.
- Variants at the edge of chromosomes where there are not enough base pairs up or downstream of the flank get shifted so that they are still the same sequence length (e.g. flank*2 + 1), but the focal base is shifted to account for the lack of base pairs either up or downstream.

## Usage

```
searchSNP --vcf output.vcf --reference reference.fasta \
          --fwd sample_R1.fastq.gz --rev sample_R2.fastq.gz \
          --outdir results --prefix sample1
```

| option | default | meaning |
|---|---|---|
| `--flank N` | 11 | bases before the variant; sequences are `2*N+1` long. Allowed 1 to 15 |
| `--threads N` | all cores | number of k-mer counting threads, split across the two read files; each file also gets one extra thread that decompresses it |
| `--outdir DIR` | `.` | output directory (created if missing) |
| `--prefix NAME` | `searchSNP` | prefix of the output files |
| `--strict` | off | stop with an error if any VCF REF disagrees with the reference FASTA |
| `--max-combos N` | 1024 | skip a variant that needs more haplotype combinations than this |
| `--ancestor-fwd FILE`, `--ancestor-rev FILE` | off | ancestral reads (recommended, see below); pairs whose ALT is found in them are left out of `<prefix>_summary_filtered.csv` |
| `--min-alt-count N` | 1 | ancestral ALT reads (forward + reverse) needed to drop a pair |

`--fwd` and `--rev` are the paired read files (typically the R1 and R2 files).

One sample per run. For many samples, loop in the shell:

```
for s in sample1 sample2; do
  searchSNP --vcf $s.vcf --reference ref.fasta --fwd ${s}_R1.fastq.gz --rev ${s}_R2.fastq.gz --prefix $s
done
```

We recommend running **searchSNP** on the ancestor reads in the same way to see if variant kmers are present in the ancestral reads.

### Ancestral filtering

A variant that was already present before the experiment is not evidence of a new mutation, even when the reads support it. Pass the ancestral reads and **searchSNP** removes the haplotype pairs whose ALT sequence is found in them:

```
searchSNP --vcf sample.vcf --reference ref.fasta \
          --fwd sample_R1.fastq.gz --rev sample_R2.fastq.gz \
          --ancestor-fwd ancestor_R1.fastq.gz --ancestor-rev ancestor_R2.fastq.gz \
          --prefix sample
```

This is optional but recommended. The ancestral reads are counted with the same sequences and k-mer index as the evolved reads, so a "pair" means the same thing in both samples and no separate ancestral run is needed.

**What is removed.** Each variant is searched as one or more haplotype pairs: the REF and ALT sequences that share the same combination of neighbouring variants (the `haplotype` column in `sequence_counts.csv`). A pair is removed **as a whole** when its ALT sequence has at least `--min-alt-count` reads in the ancestor. The REF sequence goes with it: keeping it would compare ALT reads from a haplotype where the ALT is new against REF reads from a haplotype where it was already there. Counts, proportions and diagnostics are rebuilt from the pairs that remain, and a variant with no pairs left keeps its row with blank measurement columns.

If the ancestral ALT is found but the ancestral REF of the same pair is not, the pair is still removed (the ALT is there), but treat it with more caution: that can happen when the haplotype is rare in the ancestor, and it also happens when the ALT k-mer matches somewhere else in the genome (`shared_in_reference = alt`).

With the default `--min-alt-count 1`, a single ancestral read is enough to drop a pair. A base-call error at the focal position produces exactly the ALT k-mer, so if variants are being dropped on one or two reads, raise the threshold (e.g. `--min-alt-count 3`). The `ancestor_count_fwd`/`ancestor_count_rev` columns show how much support each decision had.

Filtering answers "is this a new mutation", so it deliberately discards selection on standing variation (an ALT that was rare in the ancestor and common in the evolved sample). The unfiltered summary keeps those counts.

## Outputs

`<prefix>_summary.csv`: one row per VCF row (per ALT allele), in VCF order.

| column | meaning |
|---|---|
| `chrom, pos, ref, alt` | Copied from the VCF, so the table can be merged back onto it |
| `ref_count_fwd, alt_count_fwd, total_fwd, prop_alt_fwd` | Read counts from the forward reads |
| `ref_count_rev, alt_count_rev, total_rev, prop_alt_rev` | Same from the reverse reads |
| `window_bp` | Length (bp) of the sequences searched for this row (0 if skipped) |
| `ref_matches` | Does the VCF REF equal the reference FASTA at `pos`? (false suggests the wrong reference is being used) |
| `multiallelic`, `n_alts_at_site` | More than one ALT at this chrom/pos (from one `A,B` row or several rows) |
| `n_search_seqs` | Number of searched haplotype sequences |
| `search_status` | `full`, `shifted`, `truncated`, `skipped_too_short` or `skipped_too_many_combos` (see below) |
| `shared_kmers` | A searched sequence is also used by another position: counts may be inflated |
| `n_ambiguous_seqs` | Number of sequences excluded from the counts because they occur under both ref and alt |
| `shared_in_reference` | Which search sequences also occur OUTSIDE their focal window in the reference genome (e.g. a duplicated region), where reads from the other copy are counted too and can inflate the counts: `none`, `ref` (only a REF sequence occurs elsewhere), `alt` (only an ALT sequence does) or `both`. Empty when the check was skipped with `--no-ref-check`, or when the variant was skipped so nothing was searched. Per-sequence positions are in the `ref_occurrences*` columns below |

A sequence that is never found counts 0, and `prop_alt` is 0 when `total` is 0.

`<prefix>_summary_filtered.csv` is written only when `--ancestor-fwd`/`--ancestor-rev` are given. It has the same columns plus `ancestor_status`:

| column | meaning |
|---|---|
| `ancestor_status` | `kept` (no pair removed), `partial` (some pairs removed, counts rebuilt from the rest), `all_dropped` (every pair removed: the eight count columns, `n_search_seqs`, `shared_kmers`, `n_ambiguous_seqs` and `shared_in_reference` are blank) or empty when the variant was skipped, so there was nothing to assess |

`<prefix>_summary.csv` is still written unchanged, so the two files can be compared row for row.

`<prefix>_sequence_counts.csv` is the raw output file before summarising at each variant. This returns one row per searched sequence.

| column | meaning |
|---|---|
| `chrom, pos, ref, alt` | The focal VCF row the sequence came from |
| `type` | `ref` or `alt`, decided only by the focal variant |
| `sequence` | The forward-strand sequence searched |
| `count_fwd, count_rev` | Forward / reverse reads containing the sequence on either strand (canonical k-mer); 0 = searched but not seen |
| `window_bp, search_status, shared_kmers` | As in the summary |
| `ambiguous` | `true` if this sequence also occurs under the other label of the same variant (e.g. a long nearby insertion fills the window, so the focal variant cannot change the string). Its reads cannot tell ref from alt, so they are excluded from the summary counts |
| `ref_occurrences` | How many times this sequence occurs in the whole reference genome (any chromosome, either strand; the focal window counts as one). 0 is normal for alt sequences — the allele is not in the reference. Empty when the check was skipped or the sequence is not searchable (contains N) |
| `ref_occurrences_elsewhere` | How many of those occurrences are OUTSIDE the focal window |
| `ref_elsewhere_positions` | Where those elsewhere-occurrences are (`chrom:pos`, semicolon-separated, capped at 5) |
| `haplotype` | 1-based id of the haplotype background within the variant. The REF and ALT rows that share it are one pair, which is what ancestral filtering removes together. Restarts at 1 for each variant |
| `ancestor_count_fwd, ancestor_count_rev` | Reads in the ancestral forward / reverse files containing the sequence; 0 = searched but not seen. Empty when no ancestral reads were given, or when the sequence is not searchable |

A variant skipped by a rule (`skipped_*`) has one placeholder row with empty sequence and
empty counts, so every VCF row appears. Rows that share a k-mer show the same count.

## How the sequences are built

Every sequence is exactly `2*flank + 1` bases (23 by default):

```
[flank reference bases before POS] + [ALT allele] + [reference bases after REF, to fill the length]
```

* **SNP**: `flank` bases, the base, `flank` bases.
* **Deletion**: ALT is just the anchor base, so reference after the deleted bases is pulled in.
* **Insertion**: anchor + inserted bases, then reference. An insertion longer than `flank`
  fills the window, so only the **left junction** is searched.
* The `ref` sequence is the same window with the REF allele (the plain reference).

**Nearby variants.** Other variants whose POS is within `flank` bases are combined: every
ref/alt combination of them is searched. A sequence is labelled `ref` or `alt` **only by the
focal variant**. Other variants at the focal position are ignored, and two overlapping
variants are never put on the same haplotype.

**Chromosome ends.** Edge variants are always searched, and `search_status` says how:

| status | meaning |
|---|---|
| `full` | Normal window, `flank` bases either side |
| `shifted` | Too close to a chromosome end to centre the window, so the window slides along the chromosome and stays full length (the variant is off-centre; the sequence is still exact and just as specific). Variants near the end of the window are not combined with neighbours beyond the original `flank` range |
| `truncated` | The whole chromosome is shorter than the full window, so a shorter window (the chromosome length, at least 15 bp) is used. Shorter sequences are less specific, so treat these rows with caution |
| `skipped_too_short` | The window would be under 15 bp (or the full window, if `2*flank+1` is already below 15); nothing searched, counts are 0 |
| `skipped_too_many_combos` | More haplotype combinations than `--max-combos`; nothing searched, counts are 0 |

Filter on `search_status`: a row with 0 counts and a `skipped_*` status means "not searched", not "no reads".

**Multi-allelic rows** (`ALT = G,T`) become one row per ALT. The REF count is shared, each
ALT is counted on its own, and `multiallelic` is set so you can collapse them afterwards.

## Limitations

* **Exact matches only**, to keep false positives low. One sequencing error inside the window means that read is not counted, so absolute counts may be slightly low. Both ref and alt are affected alike, so the proportion is little changed, but low counts will be noisy.
* A read counts only if it covers the variant with `flank` bases on both sides.
* Sequences that also occur elsewhere in the genome (repeats, rRNA operons) inflate counts.
* Counts depend on the reads passing the exact-match test, so they will not always equal the variant caller's AD/AF.
* Only A/C/G/T are searched; sequences containing N are skipped. Symbolic ALTs (`<DEL>`, `*`) are skipped with a warning.
* FASTQ files must have standard 4-line records.

## Tests

```
cargo test
```

Unit tests cover sequence building (SNPs, deletions, insertions, clusters, chromosome ends) and k-mer counting; an end-to-end test simulates reads with known ref/alt counts.

## Validation

**searchSNP** has been validated on a [_Pseudomonas fluorescens_ AB1](https://www.ncbi.nlm.nih.gov/nuccore/CP148752.1) reference genome (6,759,727 bp, single contig). We created VCFs of SNPs and indels, built the corresponding genomes with **bcftools consensus**, simulated reads with **art modern**, and ran **searchSNP** using the original _Pseudomonas_ genome, the simulated VCF and the simulated reads.

The simulated VCFs include:

- Haplotype 1 with seven variants: pos=8 T>A, pos=650000 G>A, pos=1500000 T>C, pos=3200000 G>C, pos=5000000 G>C, pos=5400000 A>C, pos=6759720 A>G
- Haplotype 2 with five clustered SNPs. pos=111100 T>C, pos=111106 GA>G, pos=111110 A>G, pos=111113 A>C, pos=111124 T>C
- Haplotype 3 with 5 variants, including deletions and insertions: pos=900000 TCTTCGATCAA>T, pos=2200000 C>CGT, pos=111103 A>AAAGCCCAATAAACCACTCTGACTGGCCGAATAGGGATATAGGCAACGACATGTGCGGCGACCCTTGCGACAGTGACGCTTTCGCCGTTGCCTAAACCTAT, pos=4000000 GTAG>G, pos=5000000 G>A
- Haplotype 4 is a control with no changes.
- A phantom vcf with 1 false positive not in any of the altered genomes: pos=111113 A>T
- A mix of all haplotypes together (not including the phantom vcf).