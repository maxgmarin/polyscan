# polyscan: a tool for identifying homopolymer tracts (or near homopolymer tracts) within DNA sequences

`polyscan` is a simple Rust tool that scans a FASTA file for windows dominated by a single nucleotide. 
The user can define the window size that is scanned across the genome, the nucleotide that is searched for, and the minimum percentage required to report a window.


## Motivation

I wrote this simple tool to identify polyA tracts within gene sequences that could result in internal priming during reverse transcription of mRNA sequences for long-read RNA-sequencing.

A QC and filtering step taken by long-read isoform analysis softwares (SQUANTI, TALON) are to flag reads that start or end at polyA tracts as likely internal priming artifacts.


`polyscan` can quickly quickly process any input genome sequence and report windows of the genome which are mostly made up of a single nucleotide. This has applications for finding polyA tracts, but I set up the program to be generalizable for any nucletoide. 

`polyscan` will output a 6 column BED file which can then be used in downstream filtering and QC of long-read RNA alignments.

Paper's that dicussion this issue are linked below:
- [Nam-2002-PNAS](https://www.pnas.org/doi/10.1073/pnas.092140899)
- [Svoboda-2022-NAR-Genomics-and-Bioinformatics](https://www.pnas.org/doi/10.1073/pnas.092140899](https://academic.oup.com/nargab/article/4/2/lqac035/6592171))



## Installation

```bash
cargo install --git https://github.com/maxgmarin/polyscan
```

*(Ensure you have the [Rust toolchain](https://www.rust-lang.org/tools/install) installed.)*

## Usage

```bash
polyscan \
  --fasta <path/to/file.fasta> \
  --window-size 10 \
  --percentage 80 \
  --nucleotide A
```

To merge passing windows into tracts and write a compressed, indexed BED:

```bash
polyscan -f GRCh38.fa.gz -w 10 -p 80 -n A --mode tracts \
  -o GRCh38.polyA_tracts_W10_P80.bed.gz --index --stats \
  --bigbed GRCh38.polyA_tracts_W10_P80.bb
```

- **--fasta** / **-f**: Path to the input FASTA, plain or compressed (gzip, bgzip, bzip2, xz or zstd, detected automatically).
- **--window-size** / **-w**: Length of the sliding window (default 10).
- **--percentage** / **-p**: Minimum percentage, a whole number from 50 to 100 (default 80). A window passes when the base count is at least `ceil(p × w / 100)`. At exactly 50 with an even window size, a window can pass on both strands.
- **--nucleotide** / **-n**: Single base to detect (A, C, G, T, N). Its complement is automatically checked for the minus strand. N is its own complement, so N windows are reported on `+` only.
- **--mode** / **-m**: `windows` (default) writes one line per passing window. `tracts` merges overlapping and book-ended passing windows on the same strand into one line per tract, as `bedtools merge -s` would.
- **--output** / **-o**: Output file (default stdout). A path ending in `.gz` or `.bgz` is written as BGZF, so `tabix`, `bedtools`, IGV and pysam can read it directly.
- **--index**: Also write a tabix index (`<output>.tbi`). Needs a `.gz`/`.bgz` output.
- **--extended**: Tracts mode only. Appends two columns, described below.
- **--bigbed** *PATH*: Tracts mode only. Also writes the tracts as [bigBed](https://genome.ucsc.edu/goldenPath/help/bigBed.html), UCSC's indexed binary BED, readable by the UCSC browser, IGV, pyBigWig and bigtools. Columns match the BED output, including `--extended`. If no tracts are found, no bigBed is written and a warning is printed. The tracts are held in memory until the scan ends (about 150 MB for GRCh38 at W10/P80).
- **--threads**: Threads used to compress the bigBed (default 4).
- **--stats**: Print counts to stderr: contigs, contigs shorter than the window, passing windows, tracts and total tract length per strand, and run time.

Lowercase (soft-masked) bases are counted. Other IUPAC codes (R, Y, …) count towards no base. Contigs shorter than the window are skipped.

## Output

Both modes write **6-column BED**, sorted by contig (in FASTA order), then start, end and strand (`+` first):
1. **chrom** (contig ID)
2. **start** (0-based inclusive)
3. **end** (0-based exclusive)
4. **name** (the nucleotide being searched for)
5. **score**: in windows mode, the percentage of the target nucleotide in the window, rounded up. In tracts mode, the highest window percentage in the tract, so a stricter cutoff can be applied later without rescanning.
6. **strand**: `+` if the window or tract is a run of the target nucleotide, `-` if it is a run of its complement (the target on the reverse strand).

With `--extended`, tracts mode adds:

7. **n_windows**: number of passing windows merged into the tract
8. **span_frac**: fraction of the target nucleotide (its complement on `-`) over the whole tract, to three decimals

Example line A:
```
contig1    100     110     A   90  -
```
means `[100..110)` has 90% of `A` on the minus strand of the contig1 DNA sequence.

Example line B:
```
contig1    200     210     A   80  +
```
means `[200-210)` has 80% of `A` on the PLUS strand of the contig1 DNA sequence.

Example tract:
```
contig1    200     231     A   100 +
```
means the passing windows on the plus strand between 200 and 231 merge into one tract, and the best of them is 100% `A`.

## License & Contributing

This code is licensed under [MIT](LICENSE).

