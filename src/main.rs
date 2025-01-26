use std::fs::File;
use std::io::{BufReader, Write};
use std::error::Error;

use clap::Parser;
use bio::io::fasta;
use bio::io::bed::{Writer, Record as BedRecord};

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
    // Parse CLI
    let args = Args::parse();

    let fasta_path = args.fasta;
    let w = args.window_size;
    let p = args.percentage;
    let user_base = args.nucleotide.to_uppercase();

    // Validate user_base is exactly one char from {A,C,G,T,N}
    if user_base.len() != 1 {
        eprintln!("Error: --nucleotide must be a single character (A, C, G, T, or N).");
        std::process::exit(1);
    }
    let base_char = user_base.chars().next().unwrap();
    match base_char {
        'A' | 'C' | 'G' | 'T' | 'N' => (),
        _ => {
            eprintln!("Error: --nucleotide must be one of A, C, G, T, or N.");
            std::process::exit(1);
        }
    }

    // Validate percentage
    if p < 50.0 || p > 100.0 {
        eprintln!("Error: --percentage must be between 50.0 and 100.0");
        std::process::exit(1);
    }

    // The minimum count needed in a window to be considered "passing"
    let threshold_count: usize = ((p / 100.0) * (w as f64)).ceil() as usize;

    // We track frequencies in [A,C,G,T,N] => [0..4]
    // Let's define mapping:
    fn nuc_to_index(nuc: u8) -> Option<usize> {
        match nuc {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            b'N' | b'n' => Some(4),
            _ => None,
        }
    }


    // Complement function for one char
    fn complement_char(c: char) -> char {
        match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => 'N', // fallback
        }
    }

    // We'll find the array indices for user base + complement
    let comp_char = complement_char(base_char);
    let user_idx = nuc_to_index(base_char as u8).unwrap();
    let comp_idx = nuc_to_index(comp_char as u8).unwrap();

    // Prepare a BED writer to stdout
    let stdout = std::io::stdout();
    let mut bed_writer = Writer::new(stdout.lock());

    /// Write a BED record, placing the "strand" in aux[2].
    ///
    ///  columns: chrom, start, end, name, score, strand
    ///
    ///  - name => user base
    ///  - score => integer (rounded up) percentage
    ///  - strand => plus or minus
    fn write_bed_record<W: Write>(
        writer: &mut Writer<W>,
        chrom: &str,
        start: u64,
        end: u64,
        name: char,        // user-chosen base (not the complement)
        score_percentage: f64,  // we will round up
        strand_symbol: &str,    // e.g. "+"
    ) -> Result<(), Box<dyn Error>> {
        let mut record = BedRecord::new();

        // columns 1..3
        record.set_chrom(chrom);
        record.set_start(start);
        record.set_end(end);

        // aux[0] => name
        record.set_name(&name.to_string());

        // aux[1] => score (round up to integer)
        let ceil_int = score_percentage.ceil() as u64;
        record.set_score(&ceil_int.to_string());

        // aux[2] => strand
        record.push_aux(strand_symbol);

        writer.write(&record)?;
        Ok(())
    }

    // Open FASTA
    // let reader = fasta::Reader::new(BufReader::new(File::open(fasta_path)?));


    // Use Niffler to automatically detect compression
    let file = File::open(&fasta_path)?;
    // niffler::get_reader takes a "Box<dyn Read>", returns (reader, format)
    let (niffler_reader, _compression_format) = niffler::get_reader(Box::new(file))?;
    
    // Wrap the decompressed reader in a BufReader
    let buf = BufReader::new(niffler_reader);

    // Now create a Rust-Bio FASTA reader from that
    let reader = fasta::Reader::new(buf);





    // For each contig
    for result_record in reader.records() {
        let record = result_record?;
        let contig_id = record.id();
        let seq = record.seq();

        // Skip if contig too short
        if seq.len() < w {
            continue;
        }

        // freq array for [A,C,G,T,N]
        let mut freq = [0_usize; 5];

        // Initialize freq in the first window
        for &nuc in &seq[0..w] {
            if let Some(i) = nuc_to_index(nuc) {
                freq[i] += 1;
            }
        }

        // Check the first window
        {
            let user_count = freq[user_idx];
            let comp_count = freq[comp_idx];

            // If user base >= threshold => output plus
            if user_count >= threshold_count {
                let perc = (user_count as f64 / w as f64) * 100.0;
                write_bed_record(&mut bed_writer, contig_id, 0, w as u64, base_char, perc, "+")?;
            }

            // If complement base >= threshold => minus
            if comp_count >= threshold_count {
                let perc = (comp_count as f64 / w as f64) * 100.0;
                // We STILL label the record with the user's base, but mark strand="-"
                write_bed_record(&mut bed_writer, contig_id, 0, w as u64, base_char, perc, "-")?;
            }
        }

        // Slide the window
        for start in 1..=(seq.len() - w) {
            let leaving = seq[start - 1];
            if let Some(i) = nuc_to_index(leaving) {
                freq[i] = freq[i].saturating_sub(1);
            }

            let entering = seq[start + w - 1];
            if let Some(i) = nuc_to_index(entering) {
                freq[i] += 1;
            }

            let user_count = freq[user_idx];
            let comp_count = freq[comp_idx];

            if user_count >= threshold_count {
                let perc = (user_count as f64 / w as f64) * 100.0;
                let end = start + w;
                write_bed_record(
                    &mut bed_writer,
                    contig_id,
                    start as u64,
                    end as u64,
                    base_char,
                    perc,
                    "+"
                )?;
            }

            if comp_count >= threshold_count {
                let perc = (comp_count as f64 / w as f64) * 100.0;
                let end = start + w;
                write_bed_record(
                    &mut bed_writer,
                    contig_id,
                    start as u64,
                    end as u64,
                    base_char,
                    perc,
                    "-"
                )?;
            }
        }
    }

    Ok(())
}
