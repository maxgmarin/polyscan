//! BED output.

use std::io::{self, Write};

/// Write one BED6 line: chrom, start, end, name, score, strand.
pub fn write_window<W: Write>(
    out: &mut W,
    chrom: &str,
    start: usize,
    end: usize,
    name: char,
    score: u64,
    strand: char,
) -> io::Result<()> {
    writeln!(out, "{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}")
}
