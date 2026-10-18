//! BED output: plain text (stdout or file) or BGZF, with an optional tabix
//! index built in-process as records are written.

use std::error::Error;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::header;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_tabix as tabix;

/// True if `path` should be written as BGZF (ends in `.gz` or `.bgz`).
pub fn is_bgzf_path(path: &Path) -> bool {
    matches!(path.extension().and_then(|e| e.to_str()), Some("gz" | "bgz"))
}

/// Optional columns appended after BED6 in tracts mode with `--extended`.
pub struct Extended {
    pub n_windows: u32,
    pub span_frac: f64,
}

enum Inner {
    Plain(BufWriter<Box<dyn Write>>),
    Bgzf(bgzf::io::Writer<File>),
}

/// Destination for BED records, written in order.
pub struct Sink {
    inner: Inner,
    /// Tabix indexer and the `.tbi` path to write on `finish`.
    index: Option<(tabix::index::Indexer, PathBuf)>,
    line: String,
}

impl Sink {
    pub fn stdout() -> Self {
        let out: Box<dyn Write> = Box::new(io::stdout().lock());
        Sink::plain(out)
    }

    /// Open `path`, as BGZF if it ends in `.gz`/`.bgz`, else as plain text.
    /// With `index`, also write a tabix index to `<path>.tbi` (BGZF only).
    pub fn create(path: &Path, index: bool) -> Result<Self, Box<dyn Error>> {
        let file = File::create(path).map_err(|e| format!("cannot create {}: {e}", path.display()))?;
        if !is_bgzf_path(path) {
            return Ok(Sink::plain(Box::new(file)));
        }
        let index = index.then(|| {
            let mut indexer = tabix::index::Indexer::default();
            indexer.set_header(header::Builder::bed().build());
            let mut tbi = path.as_os_str().to_owned();
            tbi.push(".tbi");
            (indexer, PathBuf::from(tbi))
        });
        Ok(Sink { inner: Inner::Bgzf(bgzf::io::Writer::new(file)), index, line: String::new() })
    }

    fn plain(out: Box<dyn Write>) -> Self {
        Sink { inner: Inner::Plain(BufWriter::new(out)), index: None, line: String::new() }
    }

    /// Write one BED6 line (chrom, start, end, name, score, strand), plus
    /// `n_windows` and `span_frac` when `extended` is given.
    #[allow(clippy::too_many_arguments)]
    pub fn write_record(
        &mut self,
        chrom: &str,
        start: usize,
        end: usize,
        name: char,
        score: u64,
        strand: char,
        extended: Option<&Extended>,
    ) -> Result<(), Box<dyn Error>> {
        self.line.clear();
        write!(self.line, "{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")?;
        if let Some(e) = extended {
            write!(self.line, "\t{}\t{:.3}", e.n_windows, e.span_frac)?;
        }
        self.line.push('\n');

        match &mut self.inner {
            Inner::Plain(w) => w.write_all(self.line.as_bytes())?,
            Inner::Bgzf(w) => {
                let chunk_start = w.virtual_position();
                w.write_all(self.line.as_bytes())?;
                if let Some((indexer, _)) = &mut self.index {
                    let chunk = Chunk::new(chunk_start, w.virtual_position());
                    // noodles positions are 1-based inclusive: BED [start, end) is [start+1, end].
                    let first = Position::try_from(start + 1)?;
                    let last = Position::try_from(end)?;
                    indexer.add_record(chrom, first, last, chunk)?;
                }
            }
        }
        Ok(())
    }

    /// Flush the output, close the BGZF stream (writing its EOF block) and
    /// write the index if one was requested.
    pub fn finish(self) -> Result<(), Box<dyn Error>> {
        match self.inner {
            Inner::Plain(mut w) => w.flush()?,
            Inner::Bgzf(w) => {
                w.finish()?;
            }
        }
        if let Some((indexer, path)) = self.index {
            tabix::fs::write(path, &indexer.build())?;
        }
        Ok(())
    }
}
