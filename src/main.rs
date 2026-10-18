mod output;
mod scan;

use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

use bio::io::fasta;
use clap::Parser;

use output::write_window;
use scan::{complement, nuc_to_index, scan_contig, score, threshold_count};

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
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = Args::parse();

    let w = args.window_size;
    let p = args.percentage;
    let user_base = args.nucleotide.to_uppercase();

    // Validate user_base is exactly one char from {A,C,G,T,N}
    if user_base.len() != 1 {
        eprintln!("Error: --nucleotide must be a single character (A, C, G, T, or N).");
        std::process::exit(1);
    }
    let base = user_base.as_bytes()[0];
    if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
        eprintln!("Error: --nucleotide must be one of A, C, G, T, or N.");
        std::process::exit(1);
    }

    if w == 0 {
        eprintln!("Error: --window-size must be at least 1.");
        std::process::exit(1);
    }

    let threshold = threshold_count(p, w);
    let user_idx = nuc_to_index(base).unwrap();
    let comp_idx = nuc_to_index(complement(base)).unwrap();
    // N is self-complementary: checking the '-' strand would duplicate every '+' line.
    let check_comp = comp_idx != user_idx;
    let name = base as char;

    let reader = fasta::Reader::new(BufReader::new(open_input(&args.fasta)?));

    let stdout = std::io::stdout();
    let mut out = BufWriter::new(stdout.lock());

    for result_record in reader.records() {
        let record = result_record?;
        let chrom = record.id();

        for hit in scan_contig(record.seq(), w, user_idx, comp_idx) {
            let end = hit.start + w;
            // The name column is always the user's base; strand says which base ran.
            if hit.user_count >= threshold {
                write_window(&mut out, chrom, hit.start, end, name, score(hit.user_count, w), '+')?;
            }
            if check_comp && hit.comp_count >= threshold {
                write_window(&mut out, chrom, hit.start, end, name, score(hit.comp_count, w), '-')?;
            }
        }
    }

    out.flush()?;
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
    match niffler::get_reader(Box::new(File::open(path)?)) {
        Ok((reader, _format)) => Ok(reader),
        Err(niffler::Error::FileTooShort) => Ok(Box::new(File::open(path)?)),
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
