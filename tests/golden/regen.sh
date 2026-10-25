#!/usr/bin/env bash
# Regenerate the golden files listed in tests/golden/cases.tsv.
#
# windows goldens are polyscan's own output. tracts goldens are that windows
# output merged by bedtools (`merge -s -c 5,6 -o max,distinct`) and sorted the
# way polyscan sorts tracts: contigs in FASTA order, then start, end, strand.
# They are the independent reference that tracts mode must match, and are
# committed so tests don't need bedtools.
#
# Only run this when an output change is intended, and review the golden diff
# before committing.
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
command -v bedtools > /dev/null || { echo "bedtools not found" >&2; exit 1; }

grep -v '^#' tests/golden/cases.tsv | while IFS=$'\t' read -r fixture w p n modes; do
  stem="tests/golden/${fixture%.gz}.w${w}.p${p}.${n}"
  windows=$("$POLYSCAN" -f "tests/data/$fixture" -w "$w" -p "$p" -n "$n" --mode windows)
  if [[ ",$modes," == *",windows,"* ]]; then
    printf '%s' "$windows${windows:+$'\n'}" > "$stem.windows.bed"
    echo "$stem.windows.bed ($(wc -l < "$stem.windows.bed" | tr -d ' ') lines)"
  fi
  if [[ ",$modes," == *",tracts,"* && -z "$windows" ]]; then
    : > "$stem.tracts.bed"   # bedtools errors on empty input
    echo "$stem.tracts.bed (0 lines)"
  elif [[ ",$modes," == *",tracts,"* ]]; then
    printf '%s' "$windows${windows:+$'\n'}" \
      | bedtools merge -s -c 5,6 -o max,distinct -i stdin \
      | awk -v OFS='\t' -v n="$n" '!($1 in rank) { rank[$1] = ++k } { print rank[$1], $1, $2, $3, n, $4, $5 }' \
      | LC_ALL=C sort -k1,1n -k3,3n -k4,4n -k7,7 \
      | cut -f2- > "$stem.tracts.bed"
    echo "$stem.tracts.bed ($(wc -l < "$stem.tracts.bed" | tr -d ' ') lines)"
  fi
done
