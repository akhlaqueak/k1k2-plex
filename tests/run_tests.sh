#!/usr/bin/env bash
# Correctness tests for k1k2-plex. Run from the repository root (or via `make test`).
#
# For every configuration in tests/expected_counts.txt:
#   1. Ours (serial/kplex) returns the same result set as plain Bron-Kerbosch without pruning (serial/kplex-v0),
#   2. the parallel version with 1, 4 and 8 threads returns the same result set (skipped if not built),
#   3. tests/check_plexes.py independently confirms every set is a maximal (k1,k2)-plex with >= q vertices
#      and that no set is reported twice,
#   4. the number of plexes matches the expected count.
set -u
cd "$(dirname "$0")/.."
OUT=tests/out
mkdir -p "$OUT"
[ -x serial/kplex ] || { echo "serial/kplex not found: run 'make' first"; exit 1; }
[ -x serial/kplex-v0 ] || make -C serial versions > /dev/null || exit 1
PAR=parallel/kplex
[ -x "$PAR" ] || { echo "note: $PAR not built, skipping the parallel checks"; PAR=""; }

pass=0; fail=0
run() {   # run <binary> <graph> <k1> <k2> <q> <out>; extra environment via the caller
  "$1" -g "$2" -k1 "$3" -k2 "$4" -q "$5" -o "$6" > /dev/null && LC_ALL=C sort -o "$6" "$6"
}
while read -r graph k1 k2 q expected; do
  case "$graph" in ''|'#'*) continue ;; esac
  name="$(basename "$graph" .txt)-$k1-$k2-$q"
  problems=()
  run serial/kplex-v0 "$graph" "$k1" "$k2" "$q" "$OUT/$name.ref"
  run serial/kplex "$graph" "$k1" "$k2" "$q" "$OUT/$name.ours"
  cmp -s "$OUT/$name.ref" "$OUT/$name.ours" || problems+=("serial result set differs from plain BK")
  if [ -n "$PAR" ]; then
    for t in 1 4 8; do
      OMP_NUM_THREADS=$t run "$PAR" "$graph" "$k1" "$k2" "$q" "$OUT/$name.par$t"
      cmp -s "$OUT/$name.ref" "$OUT/$name.par$t" || problems+=("parallel ($t threads) result set differs")
    done
  fi
  check=$(python3 tests/check_plexes.py "$graph" "$k1" "$k2" "$q" "$OUT/$name.ours") || problems+=("$check")
  count=$(wc -l < "$OUT/$name.ours" | tr -d ' ')
  [ "$count" = "$expected" ] || problems+=("found $count plexes, expected $expected")
  if [ ${#problems[@]} -eq 0 ]; then
    pass=$((pass + 1)); printf 'PASS  %-40s %8s plexes\n' "$name" "$count"
  else
    fail=$((fail + 1)); printf 'FAIL  %-40s %s\n' "$name" "$(IFS='; '; echo "${problems[*]}")"
  fi
done < tests/expected_counts.txt
echo "$pass passed, $fail failed"
[ $fail -eq 0 ]
