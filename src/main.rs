mod bigbed;
mod output;
mod scan;
mod tracts;

use std::error::Error;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::PathBuf;
use std::time::Instant;

use bio::io::fasta::{self, FastaRead};
use clap::{Parser, ValueEnum};

use bigbed::Retained;
use output::{is_bgzf_path, Extended, Sink};
use scan::{complement, nuc_to_index, scan_contig, score, threshold_count};
use tracts::{sort_tracts, Merger, Tract};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Mode {
    /// One line per passing window
    Windows,
    /// Passing windows merged into strand-specific tracts (as `bedtools merge -s`)
    Tracts,
}

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(name = "polyscan",
          version = env!("CARGO_PKG_VERSION"),
          author = "Maximillian Marin <maximilliangmarin@gmail.com>",
          about = "Find windows in DNA sequences that have >= threshold% of a nucleotide. Outputs 6-column BED.")]
struct Args {
    /// Path to input FASTA file
    #[arg(short, long)]
    fasta: String,

    /// Length of the sliding window
    #[arg(short = 'w', long = "window-size", default_value_t = 10)]
    window_size: usize,

    /// Whole-number percentage (50-100) of the target base required in a window.
    /// A window passes when the base count is at least ceil(p * w / 100).
    /// At exactly 50 with an even window size, a window can pass on both strands.
    #[arg(short = 'p', long = "percentage", default_value_t = 80,
          value_parser = parse_percentage)]
    percentage: u32,

    /// Base to search for (A, C, G, T or N). Its complement is reported on the
    /// '-' strand. N is its own complement, so N is reported on '+' only.
    #[arg(short = 'n', long = "nucleotide", default_value = "A")]
    nucleotide: String,

    /// Output per-window lines, or merge overlapping and book-ended passing
    /// windows on the same strand into tracts. Tract score is the highest window
    /// percentage in the tract.
    #[arg(short = 'm', long = "mode", value_enum, default_value_t = Mode::Windows)]
    mode: Mode,

    /// Output file (default: stdout). A path ending in .gz or .bgz is written
    /// as BGZF; anything else as plain text.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Write a tabix index (<output>.tbi). Requires a .gz/.bgz --output.
    #[arg(long)]
    index: bool,

    /// Tracts mode only: append n_windows (passing windows merged) and
    /// span_frac (target base fraction over the whole tract) as columns 7-8.
    #[arg(long)]
    extended: bool,

    /// Tracts mode only: also write the tracts as bigBed to this path. Its
    /// columns match the BED output, including --extended.
    #[arg(long, value_name = "PATH")]
    bigbed: Option<PathBuf>,

    /// Threads for compressing the bigBed.
    #[arg(long, default_value_t = 4)]
    threads: usize,

    /// Print counts per strand to stderr.
    #[arg(long)]
    stats: bool,
}

#[derive(Default)]
struct Stats {
    contigs: u64,
    short_contigs: u64,
    windows: [u64; 2],
    tracts: [u64; 2],
    tract_bp: [u64; 2],
}

fn fail(msg: &str) -> ! {
    eprintln!("Error: {msg}");
    std::process::exit(1);
}

fn main() {
    if let Err(e) = run() {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let args = Args::parse();

    let w = args.window_size;
    let p = args.percentage;
    let user_base = args.nucleotide.to_uppercase();

    // Validate user_base is exactly one char from {A,C,G,T,N}
    if user_base.len() != 1 {
        fail("--nucleotide must be a single character (A, C, G, T, or N).");
    }
    let base = user_base.as_bytes()[0];
    if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
        fail("--nucleotide must be one of A, C, G, T, or N.");
    }
    if w == 0 {
        fail("--window-size must be at least 1.");
    }
    if args.index && !args.output.as_deref().is_some_and(is_bgzf_path) {
        fail("--index requires --output ending in .gz or .bgz.");
    }
    if args.extended && args.mode != Mode::Tracts {
        fail("--extended requires --mode tracts.");
    }
    if args.bigbed.is_some() && args.mode != Mode::Tracts {
        fail("--bigbed requires --mode tracts.");
    }
    if args.threads == 0 {
        fail("--threads must be at least 1.");
    }

    let threshold = threshold_count(p, w);
    let user_idx = nuc_to_index(base).unwrap();
    let comp_idx = nuc_to_index(complement(base)).unwrap();
    // N is self-complementary: checking the '-' strand would duplicate every '+' line.
    let check_comp = comp_idx != user_idx;
    let name = base as char;

    let mut reader = fasta::Reader::new(BufReader::new(open_input(&args.fasta)?));
    let mut sink = match &args.output {
        Some(path) => Sink::create(path, args.index)?,
        None => Sink::stdout(),
    };
    let mut stats = Stats::default();
    let mut tracts: Vec<Tract> = Vec::new();
    let mut retained = args.bigbed.as_ref().map(|_| Retained::new(args.extended));

    // One record reused for every contig, so its buffer is allocated once at the
    // size of the largest contig. `records()` allocates a fresh one per contig,
    // and on macOS the freed buffers stayed resident (3.2 GB peak on GRCh38).
    let mut record = fasta::Record::new();
    loop {
        reader.read(&mut record)?;
        if record.is_empty() {
            break;
        }
        let chrom = record.id();
        let seq = record.seq();
        stats.contigs += 1;
        if let Some(r) = &mut retained {
            r.add_contig(chrom, seq.len())?;
        }
        if seq.len() < w {
            stats.short_contigs += 1;
            continue;
        }

        match args.mode {
            Mode::Windows => {
                for hit in scan_contig(seq, w, user_idx, comp_idx) {
                    let end = hit.start + w;
                    // The name column is always the user's base; strand says which base ran.
                    if hit.user_count >= threshold {
                        stats.windows[0] += 1;
                        sink.write_record(chrom, hit.start, end, name, score(hit.user_count, w), '+', None)?;
                    }
                    if check_comp && hit.comp_count >= threshold {
                        stats.windows[1] += 1;
                        sink.write_record(chrom, hit.start, end, name, score(hit.comp_count, w), '-', None)?;
                    }
                }
            }
            Mode::Tracts => {
                let mut plus = Merger::new(w, '+');
                let mut minus = Merger::new(w, '-');
                tracts.clear();
                for hit in scan_contig(seq, w, user_idx, comp_idx) {
                    let plus_pass = hit.user_count >= threshold;
                    let minus_pass = check_comp && hit.comp_count >= threshold;
                    stats.windows[0] += plus_pass as u64;
                    stats.windows[1] += minus_pass as u64;
                    plus.push(hit.start, hit.user_count, plus_pass, &mut tracts);
                    if check_comp {
                        minus.push(hit.start, hit.comp_count, minus_pass, &mut tracts);
                    }
                }
                plus.finish(&mut tracts);
                minus.finish(&mut tracts);
                // The two strands close tracts at different times; sort so the
                // output is ordered by start, as tabix requires.
                sort_tracts(&mut tracts);

                for t in &tracts {
                    let s = (t.strand == '-') as usize;
                    stats.tracts[s] += 1;
                    stats.tract_bp[s] += (t.end - t.start) as u64;
                    let tract_score = score(t.max_count, w);
                    let target_count = if args.extended {
                        let idx = if t.strand == '+' { user_idx } else { comp_idx };
                        seq[t.start..t.end].iter().filter(|&&b| nuc_to_index(b) == Some(idx)).count()
                    } else {
                        0
                    };
                    let extended = args.extended.then(|| Extended {
                        n_windows: t.n_windows,
                        span_frac: target_count as f64 / (t.end - t.start) as f64,
                    });
                    sink.write_record(chrom, t.start, t.end, name, tract_score, t.strand, extended.as_ref())?;
                    if let Some(r) = &mut retained {
                        r.push(t, tract_score, target_count);
                    }
                }
            }
        }
    }

    sink.finish()?;
    if let (Some(r), Some(path)) = (retained, &args.bigbed) {
        if !r.write(path, name, args.threads)? {
            eprintln!("Warning: no tracts found, so no bigBed was written to {}", path.display());
        }
    }

    if args.stats {
        eprintln!("polyscan {} stats", env!("CARGO_PKG_VERSION"));
        eprintln!("contigs\t{}", stats.contigs);
        eprintln!("contigs_shorter_than_window\t{}", stats.short_contigs);
        eprintln!("threshold_count\t{threshold} of {w}");
        for (i, strand) in ["+", "-"].iter().enumerate() {
            if i == 1 && !check_comp {
                continue;
            }
            eprintln!("windows_passing\t{strand}\t{}", stats.windows[i]);
            if args.mode == Mode::Tracts {
                eprintln!("tracts\t{strand}\t{}", stats.tracts[i]);
                eprintln!("tract_bp\t{strand}\t{}", stats.tract_bp[i]);
            }
        }
        eprintln!("elapsed_seconds\t{:.2}", started.elapsed().as_secs_f64());
    }
    Ok(())
}

/// Parse `-p`: a whole number from 50 to 100. "80.0" is accepted as 80 so
/// v0.1.0 command lines keep working; "80.5" is rejected.
fn parse_percentage(s: &str) -> Result<u32, String> {
    let v: f64 = s.parse().map_err(|_| format!("'{s}' is not a number"))?;
    if v.fract() != 0.0 {
        return Err(format!("'{s}' must be a whole number"));
    }
    if !(50.0..=100.0).contains(&v) {
        return Err(format!("'{s}' must be between 50 and 100"));
    }
    Ok(v as u32)
}

/// Open the FASTA, decompressing whatever niffler detects. niffler needs a few
/// bytes to sniff the format, so a file too short for that (including an
/// empty file) is read as plain text.
fn open_input(path: &str) -> Result<Box<dyn Read>, Box<dyn Error>> {
    let open = || File::open(path).map_err(|e| format!("cannot open {path}: {e}"));
    match niffler::get_reader(Box::new(open()?)) {
        Ok((reader, _format)) => Ok(reader),
        Err(niffler::Error::FileTooShort) => Ok(Box::new(open()?)),
        Err(e) => Err(e.into()),
    }
}

#[cfg(test)]
mod tests {
    use super::parse_percentage;

    #[test]
    fn percentage_parsing() {
        assert_eq!(parse_percentage("80"), Ok(80));
        assert_eq!(parse_percentage("80.0"), Ok(80));
        assert_eq!(parse_percentage("50"), Ok(50));
        assert_eq!(parse_percentage("100"), Ok(100));
        assert!(parse_percentage("80.5").is_err());
        assert!(parse_percentage("49").is_err());
        assert!(parse_percentage("101").is_err());
        assert!(parse_percentage("abc").is_err());
    }
}
