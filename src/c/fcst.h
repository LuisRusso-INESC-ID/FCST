/*
 * Fully Compressed Suffix Tree (Russo, Navarro, Oliveira).
 *
 * A node is its suffix-array interval [l, r]; it is a leaf iff l == r.
 * The tree is represented by
 *   - an FM-index of T$ (wavelet matrix over the BWT) giving LF and Psi,
 *   - SA / ISA samples every `sa_rate` text positions,
 *   - the delta-sampled tree: a node v is sampled iff SDep(v) = 0 mod h
 *     and v = SLink^h(v') for some node v' (h = floor(delta/2)), plus the
 *     root.  For every node v there is some k < 2h <= delta such that
 *     SLink^k(v) is sampled.  The sampled tree is stored as balanced
 *     parentheses interleaved with the leaves (bitvectors B and P), with
 *     an rmM-style structure for LCSA, and SDep/h per
 *     sampled node ([l,r] is recovered from the parentheses).
 *
 * Core identity:  SDep(LCA(v,w)) = max_{0<=k<delta} k + SDep(LCSA(Psi^k(v_l), Psi^k(w_r))),
 * and LCA(v,w) = LF(letters, LCSA) for the maximising k.
 */
#ifndef FCST_H
#define FCST_H

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef struct { size_t l, r; } fcst_node;

#define FCST_NONE_POS ((size_t)-1)
static const fcst_node FCST_NONE = { FCST_NONE_POS, FCST_NONE_POS };

typedef struct fcst fcst_t;

/* Build over text[0..m-1] (no zero bytes; a unique terminator is appended).
 * delta >= 2 is the sampling parameter, sa_rate the SA/ISA sampling step.  */
fcst_t *fcst_build(const uint8_t *text, size_t m, unsigned delta, unsigned sa_rate);

/* Fast FCST (Navarro, Russo): LCA, Parent and SLink find the node as the
 * smallest LF-mapped LCSA interval, so the parenthesis sampling (delta_par)
 * needs no string depths.  SDep is stored, Elias-delta coded, only on the
 * sparser set given by delta_sdep (union-ed into the parenthesis tree).   */
fcst_t *fcst_build_fast(const uint8_t *text, size_t m, unsigned delta_par, unsigned delta_sdep,
                        unsigned sa_rate);
void    fcst_free(fcst_t *t);

size_t  fcst_n(const fcst_t *t);                 /* m + 1 (with terminator) */
size_t  fcst_sampled(const fcst_t *t);           /* sampled-tree size       */
void    fcst_space(const fcst_t *t, FILE *out);  /* space report            */
size_t  fcst_bytes(const fcst_t *t);

static inline int    fcst_is_none(fcst_node v) { return v.l == FCST_NONE_POS; }
static inline int    fcst_is_leaf(fcst_node v) { return v.l == v.r; }
static inline size_t fcst_count(fcst_node v)   { return v.r - v.l + 1; }
static inline int    fcst_eq(fcst_node a, fcst_node b) { return a.l == b.l && a.r == b.r; }
/* v is an ancestor of w (or v == w) */
static inline int    fcst_ancestor(fcst_node v, fcst_node w) { return v.l <= w.l && w.r <= v.r; }

fcst_node fcst_root(const fcst_t *t);
size_t    fcst_sdep(const fcst_t *t, fcst_node v);
fcst_node fcst_lca(const fcst_t *t, fcst_node v, fcst_node w);
fcst_node fcst_parent(const fcst_t *t, fcst_node v);
fcst_node fcst_slink(const fcst_t *t, fcst_node v);
fcst_node fcst_slink_i(const fcst_t *t, fcst_node v, size_t i);
fcst_node fcst_child(const fcst_t *t, fcst_node v, uint8_t c);
fcst_node fcst_fchild(const fcst_t *t, fcst_node v);
fcst_node fcst_nsibling(const fcst_t *t, fcst_node v);
fcst_node fcst_wlink(const fcst_t *t, fcst_node v, uint8_t c);   /* Weiner link   */
fcst_node fcst_search(const fcst_t *t, const uint8_t *p, size_t m); /* locus of p */
int       fcst_letter(const fcst_t *t, fcst_node v, size_t i);   /* i>=1; -1 = $  */

/* compressed suffix array primitives */
size_t fcst_locate(const fcst_t *t, size_t i);   /* A[i]     */
size_t fcst_isa(const fcst_t *t, size_t p);      /* A^{-1}[p] */
size_t fcst_psi(const fcst_t *t, size_t i);
size_t fcst_lf(const fcst_t *t, size_t i);

#endif
