mod output;
mod scan;

use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};

use bio::io::fasta;
use clap::Parser;

use output::write_window;
use scan::{complement, nuc_to_index, scan_contig, score, threshold_count};

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(name = "polyscan",
          version = "0.1.0",
          author = "Maximillian Marin <maximilliangmarin@gmail.com>",
          about = "Find windows in DNA sequences that have >= threshold% of a nucleotide. Outputs 6-column BED.")]
struct Args {
    /// Path to input FASTA file
    #[arg(short, long)]
    fasta: String,

    /// Window size
    #[arg(short = 'w', long = "window-size", default_value_t = 10,
          help = "Length of the sliding window")]
    window_size: usize,

    /// Percentage threshold (e.g. 80.0 for 80%)
    #[arg(short = 'p', long = "percentage", default_value_t = 80.0,
          help = "Percentage of target nucleotide required in the window")]
    percentage: f64,

    /// Single nucleotide to check (A,C,G,T,N). Its complement is automatically handled.
    #[arg(short = 'n', long = "nucleotide", default_value = "A",
          help = "nucleotide base to search for (i.e A, C, T, or G)")]
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

    if !(50.0..=100.0).contains(&p) {
        eprintln!("Error: --percentage must be between 50.0 and 100.0");
        std::process::exit(1);
    }

    let threshold = threshold_count(p, w);
    let user_idx = nuc_to_index(base).unwrap();
    let comp_idx = nuc_to_index(complement(base)).unwrap();
    let name = base as char;

    // niffler detects the input compression
    let file = File::open(&args.fasta)?;
    let (niffler_reader, _compression_format) = niffler::get_reader(Box::new(file))?;
    let reader = fasta::Reader::new(BufReader::new(niffler_reader));

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
            if hit.comp_count >= threshold {
                write_window(&mut out, chrom, hit.start, end, name, score(hit.comp_count, w), '-')?;
            }
        }
    }

    out.flush()?;
    Ok(())
}
