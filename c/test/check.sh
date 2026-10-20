#!/bin/bash
# Builds test.c with every generator (at -O0 and -O3, with all warnings),
# checks internal consistency (jump_ce(1, 16) against 2^16 calls to next(),
# jump() and long_jump() against jump_ce()), and compares outputs and jumps
# with the independent Python model in check.py.
cd "$(dirname "$0")/.." || exit 1
OUT=$(mktemp -d)
for g in "MWC128X 1 64 96" "MWC192X 2 96 144" "MWC256X 3 128 192" "GMWC128X 1 64 96 -DGMWC" "GMWC256X 3 128 192 -DGMWC" \
         "MWC64X_32 1 32 48 -DW=32" "MWC128X_32 3 64 96 -DW=32"; do
  set -- $g
  for opt in -O0 -O3; do
    cc $opt -march=native -std=c99 -Wall -Wextra -Wno-unused-function -I. -DGENFILE="\"$1.c\"" -DK=$2 -DJE=$3 -DLE=$4 $5 test/test.c -o "$OUT/t-$1$opt" || exit 1
  done
  "$OUT/t-$1-O0" > "$OUT/$1-O0.out"; "$OUT/t-$1-O3" > "$OUT/$1.out"
  cmp -s "$OUT/$1-O0.out" "$OUT/$1.out" && echo "$1: -O0 and -O3 agree" || echo "$1: -O0/-O3 DIFFER"
  tail -3 "$OUT/$1.out" | sed "s/^/$1: /"
done
python3 test/check.py "$OUT"
