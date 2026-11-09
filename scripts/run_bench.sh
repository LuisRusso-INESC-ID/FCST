#!/usr/bin/env bash
# Run the comparison benchmark on text files and append rows to the result TSVs.
#   scripts/run_bench.sh same    FILE...   # all structures on the shared Huffman+RRR CSA
#   scripts/run_bench.sh default FILE...   # each library's default configuration
# The text name (file name without extension) becomes the first column.
# Texts must not contain zero bytes.
set -uo pipefail
cd "$(dirname "$0")/.."
mode=$1; shift
BIN=build/bench_cst
case $mode in
    same)    OUT=results/results_same_csa.tsv
             STRUCTS="sada-rrr sct3-rrr sct3-t2-rrr cn-rrr cnr-rrr fcst-rrr:4 fcst-rrr:8 fcst-rrr:32 ffcst-rrr:4:32 ffcst-rrr:6:64 ffcst-rrr:8:64" ;;
    default) OUT=results/results_default.tsv
             STRUCTS="sada sct3 sct3-t2 cn cnr fcst:4 fcst:8 fcst:32 ffcst:4:32 ffcst:6:64 ffcst:8:64" ;;
    *) echo "mode must be 'same' or 'default'"; exit 1 ;;
esac
HDR="text	structure	n	bits_per_char	build_s	build_peak_MB	V1_nodes	V2_nodes	parent_ns	sdep_ns	fchild_ns	sibling_ns	child_ns	slink_ns	lca_ns"
[ -s "$OUT" ] || echo "$HDR" > "$OUT"
mkdir -p /tmp/sdsltmp
for f in "$@"; do
    name=$(basename "$f"); name=${name%.*}
    for s in $STRUCTS; do
        echo "== $name $s" >&2
        "$BIN" "$f" "$s" | grep '^RESULT' | sed "s/^RESULT/$name/" >> "$OUT" || echo "FAILED: $name $s" >&2
    done
done
