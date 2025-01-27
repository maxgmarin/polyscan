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

- **--fasta**: Path to the input FASTA.  
- **--window-size** / **-w**: Length of the sliding window (default 10).  
- **--percentage** / **-p**: Minimum % threshold (0–100, default 80).  
- **--nucleotide** / **-n**: Single base to detect (A, C, G, T, N). Its complement is automatically checked for the minus strand.

## Output

The following **6-column BED** lines will be written to stdout:
1. **chrom** (contig ID)  
2. **start** (0-based inclusive)  
3. **end** (0-based exclusive)  
4. **name** (the nucleotide being searched for)  
5. **score** (percentage of target nucleotide within window)  
6. **strand** (`+` or `-`)

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


## License & Contributing

This code is licensed under [MIT](LICENSE).

