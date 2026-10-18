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

            // BGZF decompresses to exactly the plain-text output.
            let mut text = Vec::new();
            noodles_bgzf::io::Reader::new(fs::File::open(&gz).unwrap()).read_to_end(&mut text).unwrap();
            assert_eq!(text, expected, "{}: BGZF content", gz.display());

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

    // A plain-text path writes plain text, identical to stdout.
    let c = &cases()[0];
    let plain = dir.join("out.bed");
    ok_stdout(c.run(&["-o", plain.to_str().unwrap()]));
    assert_eq!(fs::read(&plain).unwrap(), fs::read(c.golden("windows")).unwrap());
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

#[test]
fn invalid_flag_combinations_fail() {
    let fa = data("edges.fa");
    let fa = fa.to_str().unwrap();
    let dir = scratch("flags");
    let plain = dir.join("x.bed");
    for args in [
        vec!["-f", fa, "--index"],
        vec!["-f", fa, "--index", "-o", plain.to_str().unwrap()],
        vec!["-f", fa, "--extended"],
        vec!["-f", fa, "-p", "80.5"],
        vec!["-f", fa, "-p", "49"],
        vec!["-f", fa, "-w", "0"],
        vec!["-f", fa, "-n", "R"],
    ] {
        let out = Command::new(BIN).args(&args).output().unwrap();
        assert!(!out.status.success(), "expected failure: {args:?}");
    }
    fs::remove_dir_all(&dir).unwrap();
}
