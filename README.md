# FCST: Fully-Compressed Suffix Trees

Implementations of the Fully-Compressed Suffix Tree (FCST; Russo, Navarro,
Oliveira, ACM TALG 2011) and of the Fast FCST (Navarro, Russo, DCC 2014), plus
a benchmark comparing them with other compressed suffix trees.

```
src/c/              stand-alone C11 library (no dependencies)
src/sdsl/           fcst_sdsl.hpp: header-only C++ version on any SDSL CSA
tests/              test_fcst.c: checks against a naive suffix tree
bench/              bench_c.c (C library), bench_cst.cpp (comparison), csasize.cpp
scripts/            install_deps.sh, make_test_texts.py, run_bench.sh, gen_tables.py
results/            benchmark results (TSV)
```

## C library

```c
#include "fcst.h"
fcst_t *t = fcst_build(text, m, delta, sa_rate);                    /* FCST      */
fcst_t *f = fcst_build_fast(text, m, delta_par, delta_sdep, sa_rate); /* Fast FCST */
fcst_node v = fcst_root(t);           /* nodes are suffix-array intervals {l, r} */
fcst_sdep(t, v); fcst_parent(t, v); fcst_slink(t, v); fcst_lca(t, v, w);
fcst_child(t, v, 'a'); fcst_fchild(t, v); fcst_nsibling(t, v); fcst_letter(t, v, i);
fcst_wlink(t, v, 'a'); fcst_search(t, pattern, len); fcst_locate(t, i);
```

The text must not contain zero bytes; a terminator is appended.

```
make            # build/libfcst.a, build/test_fcst, build/bench_c
make test       # 144 configurations against a naive suffix tree
```

## SDSL version

```c++
#include "fcst_sdsl.hpp"
fcsts::fcst<sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<63>>, 32, 64>> t;
t.build("text.txt", delta_par, delta_sdep, /*fast=*/true);
auto v = t.parent(t.select_leaf(5));   // SDSL CST interface
```

## Benchmark

```
scripts/install_deps.sh            # SDSL-lite and CST-CN into deps/
make bench_cst                     # build/bench_cst, build/csasize
python3 scripts/make_test_texts.py # data/headers.txt, dna-rand.txt, dna-rep.txt
scripts/run_bench.sh same    data/*.txt   # -> results/results_same_csa.tsv
scripts/run_bench.sh default data/*.txt   # -> results/results_default.tsv
python3 scripts/gen_tables.py      # -> ../paper/generated/
```

`build/bench_cst FILE STRUCT verify` writes `STRUCT.verify` with nodes, depths,
children, siblings, suffix links and LCAs. These files must be identical across
structures.

Structures: `sada`, `sct3`, `sct3-t2` (SDSL), `cn`, `cnr` (CST-CN), `fcst:D`,
`ffcst:DP:DS`; the `-rrr` variants (`sada-rrr`, ..., `fcst-rrr:D`,
`ffcst-rrr:DP:DS`) all use `csa_wt<wt_huff<rrr_vector<63>>,32,64>`.

Methodology (Gog 2011; Cánovas and Navarro 2010):
- V1: ancestors of the parents of 10,000 random leaves, used for Parent, SDep,
  FChild, NSibling and Child.
- V2: suffix-link chains from 1,000 leaves, used for SLink.
- V3: 100,000 random leaf pairs, used for LCA.

V1 and V2 are capped at 100,000 nodes.

## Limitations

- Construction keeps the SA, ISA and LCP arrays in memory: about 35–50 bytes per
  character at peak.
- The `C` library's BWT uses plain bitvectors (n⌈log σ⌉ bits).
- Level ancestor queries and TDep are not implemented.
