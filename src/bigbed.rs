//! Optional bigBed output for tracts mode, written in-process with bigtools.
//!
//! bigtools needs every contig size before it starts writing, and the sizes are
//! only known once the FASTA has been read. So tracts are kept in a compact
//! form during the scan (12 bytes each, plus 8 with --extended) and the bigBed
//! is written after it. For GRCh38 at W10/P80 that is about 150 MB. With no
//! tracts, nothing is written: bigtools cannot write an empty bigBed.

use std::collections::HashMap;
use std::error::Error;
use std::path::Path;

use bigtools::beddata::BedParserStreamingIterator;
use bigtools::{BedEntry, BigBedWrite};

use crate::tracts::Tract;

const AUTOSQL_BED6: &str = r#"table polyTract
"Homopolymer tract from polyscan"
(
string chrom;      "Contig name"
uint   chromStart; "Start (0-based)"
uint   chromEnd;   "End (exclusive)"
string name;       "Target base"
uint   score;      "Max window percentage in tract"
char[1] strand;    "+ = target base run, - = complement run"
)
"#;

const AUTOSQL_EXTENDED: &str = r#"table polyTract
"Homopolymer tract from polyscan"
(
string chrom;      "Contig name"
uint   chromStart; "Start (0-based)"
uint   chromEnd;   "End (exclusive)"
string name;       "Target base"
uint   score;      "Max window percentage in tract"
char[1] strand;    "+ = target base run, - = complement run"
uint   nWindows;   "Passing windows merged"
float  spanFrac;   "Target fraction over tract span"
)
"#;

#[derive(Clone, Copy)]
struct Packed {
    start: u32,
    end: u32,
    score: u8,
    minus: bool,
}

/// `n_windows` and the target base count over the tract, for --extended.
#[derive(Clone, Copy)]
struct PackedExt {
    n_windows: u32,
    target_count: u32,
}

struct Contig {
    name: String,
    tracts: Vec<Packed>,
    ext: Vec<PackedExt>,
}

/// Tracts retained during the scan, in output order.
pub struct Retained {
    extended: bool,
    chrom_sizes: HashMap<String, u32>,
    contigs: Vec<Contig>,
}

impl Retained {
    pub fn new(extended: bool) -> Self {
        Retained { extended, chrom_sizes: HashMap::new(), contigs: Vec::new() }
    }

    /// Start a contig. Call for every contig, in FASTA order, including ones
    /// without tracts.
    pub fn add_contig(&mut self, name: &str, len: usize) -> Result<(), Box<dyn Error>> {
        let len = u32::try_from(len).map_err(|_| format!("contig {name} is too long for bigBed ({len} bp)"))?;
        self.chrom_sizes.insert(name.to_string(), len);
        self.contigs.push(Contig { name: name.to_string(), tracts: Vec::new(), ext: Vec::new() });
        Ok(())
    }

    /// Keep one tract of the current contig. Tracts must arrive in output order.
    /// `target_count` is the target base count over the tract (used with --extended).
    pub fn push(&mut self, t: &Tract, score: u64, target_count: usize) {
        let c = self.contigs.last_mut().expect("add_contig before push");
        // Fits: add_contig checked the contig length fits in u32, and score <= 100.
        c.tracts.push(Packed { start: t.start as u32, end: t.end as u32, score: score as u8, minus: t.strand == '-' });
        if self.extended {
            c.ext.push(PackedExt { n_windows: t.n_windows, target_count: target_count as u32 });
        }
    }

    /// Write the bigBed. Columns after `end` match the BED output exactly.
    /// Returns false, writing nothing, if there are no tracts: bigtools cannot
    /// write an empty bigBed.
    pub fn write(self, path: &Path, name: char, threads: usize) -> Result<bool, Box<dyn Error>> {
        let Retained { extended, chrom_sizes, contigs } = self;
        if contigs.iter().all(|c| c.tracts.is_empty()) {
            return Ok(false);
        }
        let mut writer = BigBedWrite::create_file(path, chrom_sizes)?;
        writer.autosql = Some(if extended { AUTOSQL_EXTENDED } else { AUTOSQL_BED6 }.to_string());

        let entries = contigs.into_iter().flat_map(move |c| {
            let Contig { name: chrom, tracts, ext } = c;
            tracts.into_iter().enumerate().map(move |(i, t)| {
                let strand = if t.minus { '-' } else { '+' };
                let mut rest = format!("{name}\t{}\t{strand}", t.score);
                if extended {
                    let e = ext[i];
                    let span_frac = e.target_count as f64 / (t.end - t.start) as f64;
                    rest.push_str(&format!("\t{}\t{span_frac:.3}", e.n_windows));
                }
                (chrom.clone(), BedEntry { start: t.start, end: t.end, rest })
            })
        });
        // `true`: contigs arrive in FASTA order, not sorted by name.
        let data = BedParserStreamingIterator::wrap_infallible_iter(entries, true);
        let runtime = tokio::runtime::Builder::new_multi_thread().worker_threads(threads).build()?;
        writer.write(data, runtime).map_err(|e| format!("writing bigBed {}: {e}", path.display()))?;
        Ok(true)
    }
}
