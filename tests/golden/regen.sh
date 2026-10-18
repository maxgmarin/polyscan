#!/usr/bin/env bash
# Regenerate the windows-mode golden files from tests/golden/cases.tsv.
#
# Only run this when an output change is intended, and commit the regenerated
# goldens separately from code changes so the diff shows exactly what moved.
#
#   tests/golden/regen.sh                 # uses target/release/polyscan
#   POLYSCAN=/path/to/polyscan tests/golden/regen.sh
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/../.." && pwd)
cd "$ROOT"
if [ -z "${POLYSCAN:-}" ]; then
  cargo build --release --quiet
  POLYSCAN=target/release/polyscan
fi

grep -v '^#' tests/golden/cases.tsv | while IFS=$'\t' read -r fixture w p n; do
  out="tests/golden/${fixture%.gz}.w${w}.p${p}.${n}.windows.bed"
  "$POLYSCAN" -f "tests/data/$fixture" -w "$w" -p "$p" -n "$n" > "$out"
  echo "$out ($(wc -l < "$out" | tr -d ' ') lines)"
done
