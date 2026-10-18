//! BED output: plain text (stdout or file) or BGZF, with an optional tabix
//! index built in-process as records are written.
//!
//! Files are written under a temporary name and renamed into place only once
//! every output (BED, index, bigBed) has been written, so a run that fails
//! never leaves a complete-looking partial file behind.

use std::error::Error;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{self, BufWriter, Write};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use noodles_bgzf as bgzf;
use noodles_core::Position;
use noodles_csi::binning_index::index::header;
use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
use noodles_tabix as tabix;

/// Largest position a tabix index can address (min_shift 14, depth 5).
pub const TABIX_MAX_POSITION: usize = 1 << 29;

/// True if `path` should be written as BGZF (ends in `.gz` or `.bgz`, any case).
pub fn is_bgzf_path(path: &Path) -> bool {
    path.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case("gz") || e.eq_ignore_ascii_case("bgz"))
}

/// Optional columns appended after BED6 in tracts mode with `--extended`.
pub struct Extended {
    pub n_windows: u32,
    /// Target base count over the tract, for span_frac.
    pub target_count: usize,
}

/// Append the columns after `end`: name, score, strand, and with `extended`,
/// n_windows and span_frac. Shared by the BED and bigBed writers so their
/// columns cannot drift apart.
pub fn push_rest(
    line: &mut String,
    name: char,
    score: u64,
    strand: char,
    span: usize,
    extended: Option<&Extended>,
) {
    // Writing to a String cannot fail.
    let _ = write!(line, "{name}\t{score}\t{strand}");
    if let Some(e) = extended {
        let span_frac = e.target_count as f64 / span as f64;
        let _ = write!(line, "\t{}\t{span_frac:.3}", e.n_windows);
    }
}

/// Output files written under temporary names, renamed into place by
/// `commit`. Dropping without committing removes them.
#[derive(Default)]
pub struct Staged {
    files: Vec<(PathBuf, PathBuf)>,
    committed: bool,
}

impl Staged {
    /// Register `path` and return the temporary path to write instead.
    pub fn stage(&mut self, path: &Path) -> PathBuf {
        let mut tmp = path.as_os_str().to_owned();
        tmp.push(".tmp");
        let tmp = PathBuf::from(tmp);
        self.files.push((tmp.clone(), path.to_path_buf()));
        tmp
    }

    /// Drop a staged file that turned out not to be needed.
    pub fn discard(&mut self, path: &Path) {
        self.files.retain(|(tmp, dst)| {
            let keep = dst != path;
            if !keep {
                let _ = fs::remove_file(tmp);
            }
            keep
        });
    }

    /// Rename every staged file to its final path.
    pub fn commit(mut self) -> Result<(), Box<dyn Error>> {
        for (tmp, dst) in &self.files {
            fs::rename(tmp, dst).map_err(|e| format!("cannot rename {} to {}: {e}", tmp.display(), dst.display()))?;
        }
        self.committed = true;
        Ok(())
    }
}

impl Drop for Staged {
    fn drop(&mut self) {
        if !self.committed {
            for (tmp, _) in &self.files {
                let _ = fs::remove_file(tmp);
            }
        }
    }
}

enum Inner {
    Plain(BufWriter<Box<dyn Write>>),
    Bgzf(bgzf::io::Writer<File>),
    /// Used when no index is needed: compresses blocks on several threads, but
    /// cannot report the virtual positions an index needs.
    BgzfMt(bgzf::io::MultithreadedWriter<File>),
}

/// Destination for BED records, written in order.
pub struct Sink {
    inner: Inner,
    /// Tabix indexer and the (temporary) `.tbi` path to write on `finish`.
    index: Option<(tabix::index::Indexer, PathBuf)>,
    line: String,
}

impl Sink {
    pub fn stdout() -> Self {
        let out: Box<dyn Write> = Box::new(io::stdout().lock());
        Sink::plain(out)
    }

    /// Open `path` (staged under a temporary name), as BGZF if it ends in
    /// `.gz`/`.bgz`, else as plain text. With `index`, also write a tabix index
    /// to `<path>.tbi` (BGZF only). Without an index, BGZF is compressed on
    /// `threads` threads.
    pub fn create(path: &Path, index: bool, threads: usize, staged: &mut Staged) -> Result<Self, Box<dyn Error>> {
        let tmp = staged.stage(path);
        let file = File::create(&tmp).map_err(|e| format!("cannot create {}: {e}", tmp.display()))?;
        if !is_bgzf_path(path) {
            return Ok(Sink::plain(Box::new(file)));
        }
        if !index {
            let workers = NonZeroUsize::new(threads).unwrap_or(NonZeroUsize::MIN);
            let writer = bgzf::io::MultithreadedWriter::with_worker_count(workers, file);
            return Ok(Sink { inner: Inner::BgzfMt(writer), index: None, line: String::new() });
        }
        let mut indexer = tabix::index::Indexer::default();
        indexer.set_header(header::Builder::bed().build());
        let mut tbi = path.as_os_str().to_owned();
        tbi.push(".tbi");
        let tbi_tmp = staged.stage(Path::new(&tbi));
        Ok(Sink {
            inner: Inner::Bgzf(bgzf::io::Writer::new(file)),
            index: Some((indexer, tbi_tmp)),
            line: String::new(),
        })
    }

    fn plain(out: Box<dyn Write>) -> Self {
        Sink { inner: Inner::Plain(BufWriter::with_capacity(1 << 16, out)), index: None, line: String::new() }
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
        let _ = write!(self.line, "{chrom}\t{start}\t{end}\t");
        push_rest(&mut self.line, name, score, strand, end - start, extended);
        self.line.push('\n');

        match &mut self.inner {
            Inner::Plain(w) => w.write_all(self.line.as_bytes())?,
            Inner::BgzfMt(w) => w.write_all(self.line.as_bytes())?,
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
    /// write the index if one was requested. The files stay under their
    /// temporary names until `Staged::commit`.
    pub fn finish(self) -> Result<(), Box<dyn Error>> {
        match self.inner {
            Inner::Plain(mut w) => w.flush()?,
            Inner::Bgzf(w) => {
                w.finish()?;
            }
            Inner::BgzfMt(mut w) => {
                w.finish()?;
            }
        }
        if let Some((indexer, path)) = self.index {
            tabix::fs::write(path, &indexer.build())?;
        }
        Ok(())
    }
}
