# polyscan: a tool for identifying homopolymer tracts (or near homopolymer tracts) within DNA sequences

`polyscan` is a simple Rust tool that scans a FASTA file for windows dominated by a single nucleotide, and by default merges them into tracts.
The user can define the window size that is scanned across the genome, the nucleotide that is searched for, and the minimum percentage required to report a window.


## Motivation

I wrote this simple tool to identify polyA tracts within gene sequences that could result in internal priming during reverse transcription of mRNA sequences for long-read RNA-sequencing.

A QC and filtering step taken by long-read isoform analysis softwares (SQUANTI, TALON) are to flag reads that start or end at polyA tracts as likely internal priming artifacts.


`polyscan` can quickly quickly process any input genome sequence and report windows of the genome which are mostly made up of a single nucleotide. This has applications for finding polyA tracts, but I set up the program to be generalizable for any nucletoide. 

`polyscan` will output a 6 column BED file which can then be used in downstream filtering and QC of long-read RNA alignments.

Paper's that dicussion this issue are linked below:
- [Nam-2002-PNAS](https://www.pnas.org/doi/10.1073/pnas.092140899)
- [Svoboda-2022-NAR-Genomics-and-Bioinformatics](https://academic.oup.com/nargab/article/4/2/lqac035/6592171)



## Install

### Installing prebuilt binaries

Prebuilt binaries for Linux and macOS are attached to each [GitHub release](https://github.com/maxgmarin/polyscan/releases). Each tarball holds the `polyscan` binary and this README.

```bash
# Go to the directory where you want to place the polyscan binary
cd <directory_on_PATH>

# Download the latest release (Linux x86_64; static, no runtime dependencies)
wget https://github.com/maxgmarin/polyscan/releases/download/v0.4.0/polyscan-v0.4.0-x86_64-unknown-linux-musl.tar.gz
# or macOS (Apple Silicon):
# wget https://github.com/maxgmarin/polyscan/releases/download/v0.4.0/polyscan-v0.4.0-aarch64-apple-darwin.tar.gz

# Uncompress the .tar.gz to get the binary
tar -xzf polyscan-v0.4.0-x86_64-unknown-linux-musl.tar.gz

# Test that the binary runs
./polyscan --help
./polyscan --version

# Then add the directory with the binary to your $PATH variable
```

### Building from source

Requires a [Rust toolchain](https://rustup.rs) and a C compiler (bzip2, xz, zstd and libdeflate are compiled from source).

```bash
git clone https://github.com/maxgmarin/polyscan.git
cd polyscan
cargo build --release              # binary at target/release/polyscan
cargo install --locked --path .    # or install to ~/.cargo/bin (usually on PATH)
```

A static Linux (musl) build, e.g. for an HPC cluster, is described in the [reference](docs/REFERENCE.md#build).

## Usage

```bash
polyscan \
  --fasta <path/to/file.fasta> \
  --window-size 10 \
  --percentage 80 \
  --nucleotide A
```

To write one line per passing window instead of merged tracts:

```bash
polyscan -f GRCh38.fa.gz -w 10 -p 80 -n A --mode windows
```

To write the tracts as a compressed, indexed BED and a bigBed:

```bash
polyscan -f GRCh38.fa.gz -w 10 -p 80 -n A \
  -o GRCh38.polyA_tracts_W10_P80.bed.gz --index --stats \
  --bigbed GRCh38.polyA_tracts_W10_P80.bb
```

- **--fasta** / **-f**: Path to the input FASTA, plain or compressed (gzip, bgzip, bzip2, xz or zstd, detected automatically). Use `-` to read from stdin. bgzip (BGZF) input is decompressed on several threads, so it is the fastest compressed form.
- **--window-size** / **-w**: Length of the sliding window, from 1 to 255 (default 10).
- **--percentage** / **-p**: Minimum percentage, a whole number from 50 to 100 (default 80). A window passes when the base count is at least `ceil(p × w / 100)`. At exactly 50 with an even window size, a window can pass on both strands.
- **--nucleotide** / **-n**: Single base to detect (A, C, G, T, N). Its complement is automatically checked for the minus strand. N is its own complement, so N windows are reported on `+` only.
- **--mode** / **-m**: `tracts` (default) merges overlapping and book-ended passing windows on the same strand into one line per tract, as `bedtools merge -s` would. `windows` writes one line per passing window.
- **--output** / **-o**: Output file (default stdout). A path ending in `.gz` or `.bgz` is written as BGZF, so `tabix`, `bedtools`, IGV and pysam can read it directly.
- **--index**: Also write a tabix index (`<output>.tbi`). Needs a `.gz`/`.bgz` output. Tabix cannot address positions past 2^29 bp (536,870,912), so a longer contig stops the run with an error. For such genomes, run without `--index` and index the output with `tabix -C -p bed`, which writes a CSI index.
- **--extended**: Not available with `--mode windows`. Appends two columns, described below.
- **--bigbed** *PATH*: Not available with `--mode windows`. Also writes the tracts as [bigBed](https://genome.ucsc.edu/goldenPath/help/bigBed.html), UCSC's indexed binary BED, readable by the UCSC browser, IGV, pyBigWig and bigtools. Columns match the BED output, including `--extended`. If no tracts are found, no bigBed is written and a warning is printed. The tracts are held in memory until the scan ends (about 150 MB for GRCh38 at W10/P80).
- **--threads**: Threads for decompressing BGZF input, compressing BGZF output written without `--index`, and compressing the bigBed (default 4). Indexed BGZF output is written on one thread.
- **--stats**: Print counts to stderr: contigs, contigs shorter than the window, passing windows, tracts and total tract length per strand, and run time.

Lowercase (soft-masked) bases are counted. Other IUPAC codes (R, Y, …) count towards no base. Contigs shorter than the window are skipped. Contig names must be unique and non-empty. A duplicate or empty name stops the run with an error.

Output files are written under a `.tmp` name and renamed when the run finishes, so a failed run leaves no partial output behind. A bigBed is limited to contigs of up to 2^32 bp.

## Output

Both modes write **6-column BED** (8 columns with `--extended`), sorted by contig (in FASTA order), then start, end and strand (`+` first):
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

