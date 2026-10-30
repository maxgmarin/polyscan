mod bigbed;
mod output;
mod scan;
mod tracts;

use std::collections::HashSet;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, Read, Seek, SeekFrom};
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::time::Instant;

use bio::io::fasta::{self, FastaRead};
use clap::{Parser, ValueEnum};

use bigbed::Retained;
use output::{is_bgzf_path, Extended, Sink, Staged, TABIX_MAX_POSITION};
use scan::{complement, nuc_to_index, scan_contig, score, threshold_count, MAX_WINDOW};
use tracts::{sort_tracts, Merger, Tract};

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
enum Mode {
    /// Passing windows merged into strand-specific tracts (as `bedtools merge -s`)
    Tracts,
    /// One line per passing window
    Windows,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
#[value(rename_all = "UPPER")]
enum Base {
    A,
    C,
    G,
    T,
    N,
}

impl Base {
    fn byte(self) -> u8 {
        match self {
            Base::A => b'A',
            Base::C => b'C',
            Base::G => b'G',
            Base::T => b'T',
            Base::N => b'N',
        }
    }
}

/// Command-line arguments
#[derive(Parser, Debug)]
#[command(name = "polyscan",
          version = env!("CARGO_PKG_VERSION"),
          author = "Maximillian Marin <maximilliangmarin@gmail.com>",
          about = "Find tracts in DNA sequences where windows have >= threshold% of a nucleotide. Outputs BED.")]
struct Args {
    /// Input FASTA, plain or compressed (gzip, bgzip, bzip2, xz, zstd). Use - for stdin.
    #[arg(short, long)]
    fasta: String,

    /// Length of the sliding window (1-255)
    #[arg(short = 'w', long = "window-size", default_value_t = 10)]
    window_size: usize,

    /// Whole-number percentage (50-100) of the target base required in a window.
    /// A window passes when the base count is at least ceil(p * w / 100).
    /// At exactly 50 with an even window size, a window can pass on both strands.
    #[arg(short = 'p', long = "percentage", default_value_t = 80,
          value_parser = parse_percentage)]
    percentage: u32,

    /// Base to search for. Its complement is reported on the '-' strand. N is
    /// its own complement, so N is reported on '+' only.
    #[arg(short = 'n', long = "nucleotide", value_enum, ignore_case = true, default_value = "A")]
    nucleotide: Base,

    /// Merge overlapping and book-ended passing windows on the same strand into
    /// tracts (default), or write one line per passing window. Tract score is
    /// the highest window percentage in the tract.
    #[arg(short = 'm', long = "mode", value_enum, default_value_t = Mode::Tracts)]
    mode: Mode,

    /// Output file (default: stdout). A path ending in .gz or .bgz is written
    /// as BGZF; anything else as plain text.
    #[arg(short = 'o', long = "output")]
    output: Option<PathBuf>,

    /// Write a tabix index (<output>.tbi). Requires a .gz/.bgz --output, and
    /// contigs of at most 2^29 bp (536,870,912), the most tabix can address.
    #[arg(long)]
    index: bool,

    /// Not available with --mode windows. Appends n_windows (passing windows
    /// merged) and span_frac (target base fraction over the whole tract) as
    /// columns 7-8.
    #[arg(long)]
    extended: bool,

    /// Not available with --mode windows. Also writes the tracts as bigBed to
    /// this path. Its columns match the BED output, including --extended.
    #[arg(long, value_name = "PATH")]
    bigbed: Option<PathBuf>,

    /// Threads for decompressing BGZF input, compressing BGZF output written
    /// without --index, and compressing the bigBed.
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
        // `polyscan ... | head` closes stdout early; stop quietly, like other CLI tools.
        if e.downcast_ref::<io::Error>().is_some_and(|e| e.kind() == io::ErrorKind::BrokenPipe) {
            std::process::exit(0);
        }
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let started = Instant::now();
    let args = Args::parse();

    let w = args.window_size;
    let p = args.percentage;
    let base = args.nucleotide.byte();

    if w == 0 || w > MAX_WINDOW {
        fail(&format!("--window-size must be between 1 and {MAX_WINDOW}."));
    }
    if args.index && !args.output.as_deref().is_some_and(is_bgzf_path) {
        fail("--index requires --output ending in .gz or .bgz.");
    }
    if args.extended && args.mode != Mode::Tracts {
        fail("--extended cannot be used with --mode windows.");
    }
    if args.bigbed.is_some() && args.mode != Mode::Tracts {
        fail("--bigbed cannot be used with --mode windows.");
    }
    if args.bigbed.is_some() && args.bigbed == args.output {
        fail("--bigbed and --output must be different paths.");
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

    let mut reader = fasta::Reader::new(BufReader::new(open_input(&args.fasta, args.threads)?));

    // Declared before the sink so that on an error the sink is dropped (and its
    // files closed) before the staged temporary files are removed.
    let mut staged = Staged::default();
    let mut sink = match &args.output {
        Some(path) => Sink::create(path, args.index, args.threads, &mut staged)?,
        None => Sink::stdout(),
    };
    // Create the bigBed's temporary file now, so a bad path fails before the scan.
    let bigbed_tmp = match &args.bigbed {
        Some(path) => {
            let tmp = staged.stage(path);
            File::create(&tmp).map_err(|e| format!("cannot create {}: {e}", tmp.display()))?;
            Some(tmp)
        }
        None => None,
    };

    let mut stats = Stats::default();
    let mut tracts: Vec<Tract> = Vec::new();
    let mut retained = args.bigbed.as_ref().map(|_| Retained::new(args.extended));
    let mut seen: HashSet<String> = HashSet::new();

    // One record reused for every contig, so its buffer is allocated once at the
    // size of the largest contig. `records()` allocates a fresh one per contig,
    // and on macOS the freed buffers stayed resident (3.2 GB peak on GRCh38).
    let mut record = fasta::Record::new();
    loop {
        reader.read(&mut record).map_err(|e| format!("reading {}: {e}", args.fasta))?;
        if record.is_empty() {
            break;
        }
        let chrom = record.id();
        let seq = record.seq();

        if chrom.is_empty() {
            return Err(format!("{}: a FASTA record has an empty name", args.fasta).into());
        }
        if !seen.insert(chrom.to_string()) {
            return Err(format!("{}: duplicate contig name '{chrom}'", args.fasta).into());
        }
        if args.index && seq.len() > TABIX_MAX_POSITION {
            return Err(format!(
                "contig {chrom} is {} bp, longer than the 2^29 bp a tabix index can address. \
                 Run without --index and index the output with `tabix -C -p bed` (CSI) instead.",
                seq.len()
            )
            .into());
        }

        stats.contigs += 1;
        if let Some(r) = &mut retained {
            r.add_contig(chrom, seq.len())?;
        }
        if seq.len() < w {
            stats.short_contigs += 1;
            continue;
        }

        match args.mode {
            // The closures are forced inline: called once per passing window,
            // a separate call cost more than the scan itself.
            Mode::Windows => {
                scan_contig(seq, w, user_idx, comp_idx, threshold, #[inline(always)] |hit| {
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
                    Ok::<(), Box<dyn Error>>(())
                })?;
            }
            Mode::Tracts => {
                let mut plus = Merger::new(w, '+');
                let mut minus = Merger::new(w, '-');
                tracts.clear();
                // Only passing windows are fed to the mergers: a tract left open
                // is closed by the next passing window past its end, or by `finish`.
                scan_contig(seq, w, user_idx, comp_idx, threshold, #[inline(always)] |hit| {
                    if hit.user_count >= threshold {
                        stats.windows[0] += 1;
                        plus.push(hit.start, hit.user_count, true, &mut tracts);
                    }
                    if check_comp && hit.comp_count >= threshold {
                        stats.windows[1] += 1;
                        minus.push(hit.start, hit.comp_count, true, &mut tracts);
                    }
                    Ok::<(), Box<dyn Error>>(())
                })?;
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
                    let extended = args.extended.then_some(Extended { n_windows: t.n_windows, target_count });
                    sink.write_record(chrom, t.start, t.end, name, tract_score, t.strand, extended.as_ref())?;
                    if let Some(r) = &mut retained {
                        r.push(t, tract_score, target_count);
                    }
                }
            }
        }
    }

    sink.finish()?;
    if let (Some(r), Some(tmp), Some(path)) = (retained, &bigbed_tmp, &args.bigbed) {
        if !r.write(tmp, name, args.threads)? {
            staged.discard(path);
            eprintln!("Warning: no tracts found, so no bigBed was written to {}", path.display());
        }
    }
    // Every output is complete: move them all into place.
    staged.commit()?;

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
/// v0.1.0 command lines keep working; "80.5" and "8e1" are rejected.
fn parse_percentage(s: &str) -> Result<u32, String> {
    let (whole, frac) = s.split_once('.').unwrap_or((s, ""));
    if whole.is_empty() || !whole.bytes().all(|b| b.is_ascii_digit()) || !frac.bytes().all(|b| b == b'0') {
        return Err(format!("'{s}' must be a whole number"));
    }
    match whole.parse::<u32>() {
        Ok(v) if (50..=100).contains(&v) => Ok(v),
        _ => Err(format!("'{s}' must be between 50 and 100")),
    }
}

/// The 28-byte empty block that ends every complete BGZF file.
const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03,
    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Open the FASTA (`-` for stdin). The input is opened once and its first
/// bytes are peeked, so pipes and process substitutions work.
///
/// BGZF input (bgzip, the usual form of a compressed reference with a
/// .fai/.gzi) is decompressed on `threads` threads. Anything else goes through
/// niffler, which detects gzip, bzip2, xz and zstd and decompresses on one
/// thread. Input too short for niffler to sniff (under 5 bytes, including an
/// empty file) is read as plain text.
fn open_input(path: &str, threads: usize) -> Result<Box<dyn Read>, Box<dyn Error>> {
    let mut raw: Box<dyn Read + Send> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path).map_err(|e| format!("cannot open {path}: {e}"))?)
    };

    let mut head = Vec::with_capacity(16);
    (&mut raw).take(16).read_to_end(&mut head).map_err(|e| format!("reading {path}: {e}"))?;
    let is_bgzf = head.len() == 16 && head[..4] == [0x1f, 0x8b, 0x08, 0x04] && head[12..14] == *b"BC";
    let short = head.len() < 5;
    let input = io::Cursor::new(head).chain(raw);

    if is_bgzf {
        if path != "-" {
            warn_if_no_bgzf_eof(path);
        }
        let workers = NonZeroUsize::new(threads).unwrap_or(NonZeroUsize::MIN);
        return Ok(Box::new(noodles_bgzf::io::MultithreadedReader::with_worker_count(workers, input)));
    }
    if short {
        return Ok(Box::new(input));
    }
    let (reader, _format) = niffler::get_reader(Box::new(input)).map_err(|e| format!("reading {path}: {e}"))?;
    Ok(reader)
}

/// Warn, as htslib does, if a BGZF file lacks its EOF block: it was probably
/// truncated, and a cut at a block boundary would otherwise go unnoticed.
fn warn_if_no_bgzf_eof(path: &str) {
    let has_eof = (|| -> io::Result<bool> {
        let mut f = File::open(path)?;
        if f.metadata()?.len() < BGZF_EOF.len() as u64 {
            return Ok(false);
        }
        f.seek(SeekFrom::End(-(BGZF_EOF.len() as i64)))?;
        let mut tail = [0u8; 28];
        f.read_exact(&mut tail)?;
        Ok(tail == BGZF_EOF)
    })();
    // Only a definite "no" warns; a file we cannot seek (e.g. a FIFO) is skipped.
    if let Ok(false) = has_eof {
        eprintln!("Warning: {path} has no BGZF EOF marker; it may be truncated");
    }
}

#[cfg(test)]
mod tests {
    use super::parse_percentage;

    #[test]
    fn percentage_parsing() {
        assert_eq!(parse_percentage("80"), Ok(80));
        assert_eq!(parse_percentage("80.0"), Ok(80));
        assert_eq!(parse_percentage("80.00"), Ok(80));
        assert_eq!(parse_percentage("50"), Ok(50));
        assert_eq!(parse_percentage("100"), Ok(100));
        for bad in ["80.5", "49", "101", "abc", "8e1", "1e2", "", ".0", "-80", "+80", " 80"] {
            assert!(parse_percentage(bad).is_err(), "{bad:?} should be rejected");
        }
    }
}
