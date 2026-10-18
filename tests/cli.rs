//! Runs the polyscan binary on each case in tests/golden/cases.tsv and checks
//! stdout byte for byte against the committed golden file.
//!
//! To change a golden deliberately, run tests/golden/regen.sh and commit the
//! result separately from the code change.

use std::fs;
use std::path::Path;
use std::process::Command;

const BIN: &str = env!("CARGO_BIN_EXE_polyscan");

struct Case {
    fixture: String,
    window: String,
    percentage: String,
    nucleotide: String,
}

fn cases() -> Vec<Case> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let text = fs::read_to_string(root.join("tests/golden/cases.tsv")).unwrap();
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 4, "bad cases.tsv line: {l:?}");
            Case {
                fixture: f[0].to_string(),
                window: f[1].to_string(),
                percentage: f[2].to_string(),
                nucleotide: f[3].to_string(),
            }
        })
        .collect()
}

#[test]
fn windows_mode_matches_golden() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let cases = cases();
    assert!(!cases.is_empty());

    let mut failures = Vec::new();
    for c in &cases {
        let golden_name = format!(
            "{}.w{}.p{}.{}.windows.bed",
            c.fixture.trim_end_matches(".gz"),
            c.window,
            c.percentage,
            c.nucleotide
        );
        let expected = fs::read(root.join("tests/golden").join(&golden_name)).unwrap();

        let out = Command::new(BIN)
            .arg("-f")
            .arg(root.join("tests/data").join(&c.fixture))
            .args(["-w", &c.window, "-p", &c.percentage, "-n", &c.nucleotide])
            .output()
            .unwrap();

        if !out.status.success() {
            failures.push(format!(
                "{golden_name}: exit {:?}, stderr: {}",
                out.status.code(),
                String::from_utf8_lossy(&out.stderr)
            ));
        } else if out.stdout != expected {
            failures.push(format!("{golden_name}: stdout differs from golden"));
        }
    }
    assert!(failures.is_empty(), "\n{}", failures.join("\n"));
}
