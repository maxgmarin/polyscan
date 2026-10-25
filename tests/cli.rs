//! Runs the polyscan binary on each case in tests/golden/cases.tsv and checks
//! stdout byte for byte against the committed golden files, then checks BGZF,
//! tabix and --extended output.
//!
//! To change a golden deliberately, run tests/golden/regen.sh and review the
//! diff before committing.

use std::fs;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const BIN: &str = env!("CARGO_BIN_EXE_polyscan");

fn root() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn data(name: &str) -> PathBuf {
    root().join("tests/data").join(name)
}

struct Case {
    fixture: String,
    window: String,
    percentage: String,
    nucleotide: String,
    modes: Vec<String>,
}

impl Case {
    fn golden(&self, mode: &str) -> PathBuf {
        root().join("tests/golden").join(format!(
            "{}.w{}.p{}.{}.{mode}.bed",
            self.fixture.trim_end_matches(".gz"),
            self.window,
            self.percentage,
            self.nucleotide
        ))
    }

    fn run(&self, extra: &[&str]) -> Output {
        Command::new(BIN)
            .arg("-f")
            .arg(data(&self.fixture))
            .args(["-w", &self.window, "-p", &self.percentage, "-n", &self.nucleotide])
            .args(extra)
            .output()
            .unwrap()
    }
}

fn cases() -> Vec<Case> {
    let text = fs::read_to_string(root().join("tests/golden/cases.tsv")).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 5, "bad cases.tsv line: {l:?}");
            Case {
                fixture: f[0].to_string(),
                window: f[1].to_string(),
                percentage: f[2].to_string(),
                nucleotide: f[3].to_string(),
                modes: f[4].split(',').map(String::from).collect(),
            }
        })
        .collect()
}

fn ok_stdout(out: Output) -> Vec<u8> {
    assert!(out.status.success(), "exit {:?}: {}", out.status.code(), String::from_utf8_lossy(&out.stderr));
    out.stdout
}

/// A fresh, empty scratch directory for one test.
fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("polyscan-test-{}-{name}", std::process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn stdout_matches_golden() {
    let cases = cases();
    assert!(!cases.is_empty());

    let mut failures = Vec::new();
    for c in &cases {
        for mode in &c.modes {
            let golden = c.golden(mode);
            let expected = fs::read(&golden).unwrap();
            let out = c.run(&["--mode", mode]);
            let name = golden.file_name().unwrap().to_string_lossy().into_owned();
            if !out.status.success() {
                failures.push(format!("{name}: exit {:?}: {}", out.status.code(), String::from_utf8_lossy(&out.stderr)));
            } else if out.stdout != expected {
                failures.push(format!("{name}: stdout differs from golden"));
            }
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}

#[test]
fn bgzf_output_and_tabix_index() {
    let dir = scratch("bgzf");
    for c in cases().iter().filter(|c| c.fixture.starts_with("random_1mb") || c.fixture.starts_with("multi")) {
        for mode in &c.modes {
            let expected = fs::read(c.golden(mode)).unwrap();
            let gz = dir.join(format!("out.{mode}.bed.gz"));
            ok_stdout(c.run(&["--mode", mode, "-o", gz.to_str().unwrap(), "--index"]));

            // BGZF decompresses to exactly the plain-text output, both indexed
            // (single-threaded writer) and not (multithreaded writer).
            let mut text = Vec::new();
            noodles_bgzf::io::Reader::new(fs::File::open(&gz).unwrap()).read_to_end(&mut text).unwrap();
            assert_eq!(text, expected, "{}: BGZF content", gz.display());
            let mt = dir.join(format!("mt.{mode}.bed.bgz"));
            ok_stdout(c.run(&["--mode", mode, "-o", mt.to_str().unwrap(), "--threads", "4"]));
            let mut text = Vec::new();
            noodles_bgzf::io::Reader::new(fs::File::open(&mt).unwrap()).read_to_end(&mut text).unwrap();
            assert_eq!(text, expected, "{}: multithreaded BGZF content", mt.display());
            assert!(!dir.join(format!("mt.{mode}.bed.bgz.tbi")).exists());

            // Region queries through the index return exactly the overlapping records.
            let lines: Vec<&str> = std::str::from_utf8(&expected).unwrap().lines().collect();
            let mut reader = noodles_tabix::io::indexed_reader::Builder::default().build_from_path(&gz).unwrap();
            let chroms: Vec<&str> = lines.iter().map(|l| l.split('\t').next().unwrap()).collect();
            for (chrom, from, to) in [(chroms[0], 1, 1_000_000), (chroms[0], 5_000, 5_400), (chroms[0], 250_001, 260_000), (chroms[chroms.len() - 1], 1, 30)] {
                let region = format!("{chrom}:{from}-{to}").parse().unwrap();
                let got: Vec<String> =
                    reader.query(&region).unwrap().map(|r| r.unwrap().as_ref().to_string()).collect();
                // 1-based [from, to] overlaps BED [start, end) when start < to && end >= from.
                let want: Vec<String> = lines
                    .iter()
                    .filter(|l| {
                        let f: Vec<&str> = l.split('\t').collect();
                        let (s, e): (usize, usize) = (f[1].parse().unwrap(), f[2].parse().unwrap());
                        f[0] == chrom && s < to && e >= from
                    })
                    .map(|l| l.to_string())
                    .collect();
                assert_eq!(got, want, "{} {chrom}:{from}-{to}", gz.display());
            }
        }
    }

    // A plain-text path writes plain text, identical to stdout. With no
    // --mode, the output is tracts.
    let c = cases().into_iter().find(|c| c.fixture == "bookended.fa").unwrap();
    assert_ne!(fs::read(c.golden("tracts")).unwrap(), fs::read(c.golden("windows")).unwrap());
    let plain = dir.join("out.bed");
    ok_stdout(c.run(&["-o", plain.to_str().unwrap()]));
    assert_eq!(fs::read(&plain).unwrap(), fs::read(c.golden("tracts")).unwrap());
    fs::remove_dir_all(&dir).unwrap();
}

/// Brute-force check of the --extended columns on the random contig.
#[test]
fn extended_columns_match_recount() {
    let (rdr, _) = niffler::get_reader(Box::new(fs::File::open(data("random_1mb.fa.gz")).unwrap())).unwrap();
    let rec = bio::io::fasta::Reader::new(BufReader::new(rdr)).records().next().unwrap().unwrap();
    let seq: Vec<u8> = rec.seq().to_ascii_uppercase();

    for c in cases().iter().filter(|c| c.fixture == "random_1mb.fa.gz") {
        let w: usize = c.window.parse().unwrap();
        let p: usize = c.percentage.parse().unwrap();
        let threshold = (p * w).div_ceil(100);
        let base = c.nucleotide.as_bytes()[0];
        let comp = match base { b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C', _ => b'N' };

        let out = String::from_utf8(ok_stdout(c.run(&["--mode", "tracts", "--extended"]))).unwrap();
        let golden = fs::read_to_string(c.golden("tracts")).unwrap();
        assert_eq!(out.lines().count(), golden.lines().count());

        for (line, gold) in out.lines().zip(golden.lines()) {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), 8, "{line}");
            assert_eq!(f[..6].join("\t"), gold, "BED6 part matches plain tracts output");
            let (start, end): (usize, usize) = (f[1].parse().unwrap(), f[2].parse().unwrap());
            let target = if f[5] == "+" { base } else { comp };
            let count = |s: usize, e: usize| seq[s..e].iter().filter(|&&b| b == target).count();

            let n_windows = (start..=end - w).filter(|&s| count(s, s + w) >= threshold).count();
            assert_eq!(f[6], n_windows.to_string(), "n_windows: {line}");
            let span_frac = format!("{:.3}", count(start, end) as f64 / (end - start) as f64);
            assert_eq!(f[7], span_frac, "span_frac: {line}");
        }
    }
}

/// The bigBed read back through bigtools matches the BED output record for
/// record, with and without --extended, and region queries return the
/// overlapping subset.
#[test]
fn bigbed_matches_bed_output() {
    use bigtools::BigBedRead;

    let dir = scratch("bigbed");
    for c in cases().iter().filter(|c| c.fixture.starts_with("random_1mb") || c.fixture == "multi_contig.fa") {
        for extended in [false, true] {
            let bb = dir.join(format!("out.{}.{extended}.bb", c.window));
            let mut args = vec!["--mode", "tracts", "--bigbed", bb.to_str().unwrap()];
            if extended {
                args.push("--extended");
            }
            let bed = String::from_utf8(ok_stdout(c.run(&args))).unwrap();
            let lines: Vec<&str> = bed.lines().collect();
            assert!(!lines.is_empty());

            let mut reader = BigBedRead::open_file(&bb).unwrap();
            assert_eq!(reader.item_count().unwrap(), lines.len() as u64, "{}", bb.display());
            let autosql = reader.autosql().unwrap().unwrap();
            assert!(autosql.contains("polyTract"));
            assert_eq!(autosql.contains("spanFrac"), extended);

            // Every record, contig by contig in FASTA order.
            let mut chroms: Vec<&str> = Vec::new();
            for l in &lines {
                let chrom = l.split('\t').next().unwrap();
                if chroms.last() != Some(&chrom) {
                    chroms.push(chrom);
                }
            }
            let mut got = Vec::new();
            for chrom in &chroms {
                let len = reader.chroms().iter().find(|ci| ci.name == *chrom).unwrap().length;
                for e in reader.get_interval(chrom, 0, len).unwrap() {
                    let e = e.unwrap();
                    got.push(format!("{chrom}\t{}\t{}\t{}", e.start, e.end, e.rest));
                }
            }
            assert_eq!(got, lines, "{}: all records", bb.display());

            // A region query returns the records overlapping [from, to).
            let chrom = chroms[0];
            let (from, to) = (5_000u32, 60_000u32);
            let got: Vec<String> = reader
                .get_interval(chrom, from, to)
                .unwrap()
                .map(|e| e.unwrap())
                .map(|e| format!("{chrom}\t{}\t{}\t{}", e.start, e.end, e.rest))
                .collect();
            let want: Vec<&str> = lines
                .iter()
                .copied()
                .filter(|l| {
                    let f: Vec<&str> = l.split('\t').collect();
                    let (s, e): (u32, u32) = (f[1].parse().unwrap(), f[2].parse().unwrap());
                    f[0] == chrom && s < to && e > from
                })
                .collect();
            assert_eq!(got, want, "{}: region {chrom}:{from}-{to}", bb.display());
        }
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn bigbed_skipped_when_no_tracts() {
    let dir = scratch("bigbed-empty");
    let bb = dir.join("empty.bb");
    let out = Command::new(BIN)
        .arg("-f")
        .arg(data("short_contig.fa"))
        .args(["--mode", "tracts", "--bigbed", bb.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("no bigBed was written"));
    assert!(!bb.exists());
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn invalid_flag_combinations_fail() {
    let fa = data("edges.fa");
    let fa = fa.to_str().unwrap();
    let dir = scratch("flags");
    let plain = dir.join("x.bed");
    let bb = dir.join("x.bb");
    let (plain, bb) = (plain.to_str().unwrap(), bb.to_str().unwrap());
    for args in [
        vec!["-f", fa, "--index"],
        vec!["-f", fa, "--index", "-o", plain],
        vec!["-f", fa, "--mode", "windows", "--extended"],
        vec!["-f", fa, "-p", "80.5"],
        vec!["-f", fa, "-p", "49"],
        vec!["-f", fa, "-p", "8e1"],
        vec!["-f", fa, "-w", "0"],
        vec!["-f", fa, "-n", "R"],
        vec!["-f", fa, "-n", "AT"],
        vec!["-f", fa, "--mode", "windows", "--bigbed", bb],
        vec!["-f", fa, "--mode", "tracts", "--bigbed", bb, "--threads", "0"],
        vec!["-f", fa, "--mode", "tracts", "--bigbed", bb, "-o", bb],
    ] {
        let out = Command::new(BIN).args(&args).output().unwrap();
        assert!(!out.status.success(), "expected failure: {args:?}");
    }
    assert_eq!(fs::read_dir(&dir).unwrap().count(), 0, "a rejected command line wrote files");
    fs::remove_dir_all(&dir).unwrap();
}

/// Write `src` as BGZF (as bgzip would), optionally without the EOF block.
fn write_bgzf(src: &Path, dst: &Path, with_eof: bool) {
    let mut w = noodles_bgzf::io::Writer::new(fs::File::create(dst).unwrap());
    std::io::Write::write_all(&mut w, &fs::read(src).unwrap()).unwrap();
    w.finish().unwrap();
    if !with_eof {
        let bytes = fs::read(dst).unwrap();
        fs::write(dst, &bytes[..bytes.len() - 28]).unwrap();
    }
}

fn plain_text_of(name: &str) -> Vec<u8> {
    let (mut r, _) = niffler::get_reader(Box::new(fs::File::open(data(name)).unwrap())).unwrap();
    let mut v = Vec::new();
    r.read_to_end(&mut v).unwrap();
    v
}

/// BGZF input goes through the multithreaded reader; output must match the
/// goldens for any thread count.
#[test]
fn bgzf_input_matches_golden() {
    let dir = scratch("bgzf-in");
    for c in cases().iter().filter(|c| c.fixture == "multi_contig.fa" || c.fixture.starts_with("random_1mb")) {
        let plain = dir.join("in.fa");
        fs::write(&plain, plain_text_of(&c.fixture)).unwrap();
        let bgz = dir.join("in.fa.gz");
        write_bgzf(&plain, &bgz, true);
        let expected = fs::read(c.golden("tracts")).unwrap();
        for threads in ["1", "4"] {
            let out = Command::new(BIN)
                .arg("-f")
                .arg(&bgz)
                .args(["-w", &c.window, "-p", &c.percentage, "-n", &c.nucleotide])
                .args(["--mode", "tracts", "--threads", threads])
                .output()
                .unwrap();
            assert!(out.stderr.is_empty(), "{}", String::from_utf8_lossy(&out.stderr));
            assert_eq!(ok_stdout(out), expected, "{} threads={threads}", c.fixture);
        }
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn truncated_bgzf_input_warns() {
    let dir = scratch("bgzf-trunc");
    let bgz = dir.join("t.fa.gz");
    write_bgzf(&data("multi_contig.fa"), &bgz, false);
    let out = Command::new(BIN).arg("-f").arg(&bgz).output().unwrap();
    assert!(out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("no BGZF EOF marker"));
    fs::remove_dir_all(&dir).unwrap();
}

/// stdin (`-f -`), plain, gzip and BGZF, gives the same output as the file.
#[test]
fn stdin_input() {
    let dir = scratch("stdin");
    let plain = data("multi_contig.fa");
    let bgz = dir.join("in.fa.gz");
    write_bgzf(&plain, &bgz, true);
    let expected = fs::read(root().join("tests/golden/multi_contig.fa.w10.p80.A.tracts.bed")).unwrap();
    for input in [plain, data("multi_contig.fa.gz"), bgz] {
        let out = Command::new(BIN)
            .args(["-f", "-", "--mode", "tracts"])
            .stdin(fs::File::open(&input).unwrap())
            .output()
            .unwrap();
        assert_eq!(ok_stdout(out), expected, "stdin from {}", input.display());
    }
    fs::remove_dir_all(&dir).unwrap();
}

/// Bad contig names stop the run, and no output file (final or temporary) is left.
#[test]
fn bad_contig_names_fail_cleanly() {
    let dir = scratch("names");
    for (label, fasta) in [
        ("duplicate", ">c1\nAAAAAAAAAAAA\n>c2\nAAAAAAAAAAAA\n>c1\nAAAAAAAAAAAA\n"),
        ("empty", ">c1\nAAAAAAAAAAAA\n>\nAAAAAAAAAAAA\n"),
    ] {
        let fa = dir.join("in.fa");
        fs::write(&fa, fasta).unwrap();
        let out = Command::new(BIN)
            .arg("-f")
            .arg(&fa)
            .args(["--mode", "tracts", "--index", "-o"])
            .arg(dir.join("out.bed.gz"))
            .arg("--bigbed")
            .arg(dir.join("out.bb"))
            .output()
            .unwrap();
        assert!(!out.status.success(), "{label} name accepted");
        let left: Vec<_> = fs::read_dir(&dir).unwrap().map(|e| e.unwrap().file_name()).collect();
        assert_eq!(left, vec!["in.fa"], "{label}: files left behind");
    }
    fs::remove_dir_all(&dir).unwrap();
}

#[test]
fn lowercase_nucleotide_and_stats() {
    let c = &cases().into_iter().find(|c| c.fixture == "multi_contig.fa").unwrap();
    let upper = c.run(&["--mode", "tracts"]);
    let lower = Command::new(BIN)
        .arg("-f")
        .arg(data(&c.fixture))
        .args(["-n", "a", "--mode", "tracts", "--stats"])
        .output()
        .unwrap();
    assert!(lower.status.success());
    assert_eq!(lower.stdout, upper.stdout);

    let stats = String::from_utf8(lower.stderr).unwrap();
    let golden = fs::read_to_string(c.golden("tracts")).unwrap();
    let plus = golden.lines().filter(|l| l.ends_with("\t+")).count();
    let minus = golden.lines().filter(|l| l.ends_with("\t-")).count();
    for line in ["contigs\t4", "contigs_shorter_than_window\t1", "threshold_count\t8 of 10"] {
        assert!(stats.lines().any(|l| l == line), "missing {line:?} in:\n{stats}");
    }
    assert!(stats.lines().any(|l| l == format!("tracts\t+\t{plus}")), "{stats}");
    assert!(stats.lines().any(|l| l == format!("tracts\t-\t{minus}")), "{stats}");
}
