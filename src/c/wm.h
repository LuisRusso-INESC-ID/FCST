/* Wavelet matrix over a byte alphabet (used for the BWT / FM-index). */
#ifndef FCST_WM_H
#define FCST_WM_H

#include "bits.h"

typedef struct {
    size_t   n;
    unsigned L;     /* levels = ceil(log2 sigma) */
    bv_t    *lv;
    size_t  *z;     /* zeros per level           */
} wm_t;

void     wm_build(wm_t *w, const uint8_t *s, size_t n, unsigned sigma);
void     wm_free(wm_t *w);
unsigned wm_access(const wm_t *w, size_t i);
size_t   wm_rank(const wm_t *w, unsigned c, size_t i);            /* occ of c in [0,i) */
unsigned wm_access_rank(const wm_t *w, size_t i, size_t *rank);   /* c=s[i], rank_c(i) */
size_t   wm_select(const wm_t *w, unsigned c, size_t k);          /* k-th c, k >= 1    */
size_t   wm_bytes(const wm_t *w);

#endif
