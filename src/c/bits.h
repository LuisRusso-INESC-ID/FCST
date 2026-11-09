/* Plain bitvectors with rank/select, and packed integer arrays. */
#ifndef FCST_BITS_H
#define FCST_BITS_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    uint64_t *w;        /* bits, LSB-first inside each word            */
    size_t    n, nw;    /* length in bits, number of 64-bit words      */
    uint64_t *sb;       /* sb[j] = rank1 at bit 512*j, j = 0..nsb      */
    size_t    nsb;
    uint32_t *s1, *s0;  /* superblock holding every 4096-th one / zero */
    size_t    ns1, ns0;
    size_t    ones;
} bv_t;

void   bv_init(bv_t *b, size_t n);           /* all zeros            */
void   bv_build(bv_t *b);                    /* call after setting    */
void   bv_free(bv_t *b);
size_t bv_rank1(const bv_t *b, size_t i);    /* ones in [0, i)        */
size_t bv_select1(const bv_t *b, size_t k);  /* pos of k-th one, k>=1 */
size_t bv_select0(const bv_t *b, size_t k);  /* pos of k-th zero      */
size_t bv_bytes(const bv_t *b);

static inline size_t bv_rank0(const bv_t *b, size_t i) { return i - bv_rank1(b, i); }
static inline void   bv_set(bv_t *b, size_t i) { b->w[i >> 6] |= 1ULL << (i & 63); }
static inline int    bv_get(const bv_t *b, size_t i) { return (int)((b->w[i >> 6] >> (i & 63)) & 1); }

typedef struct {
    uint64_t *w;
    size_t    n;
    unsigned  width;
} pa_t;

void     pa_init(pa_t *a, size_t n, unsigned width);
void     pa_free(pa_t *a);
void     pa_set(pa_t *a, size_t i, uint64_t v);
size_t   pa_bytes(const pa_t *a);
unsigned bits_needed(uint64_t maxval);      /* >= 1 */

static inline uint64_t pa_get(const pa_t *a, size_t i)
{
    size_t   pos = i * a->width, q = pos >> 6;
    unsigned off = pos & 63;
    uint64_t v   = a->w[q] >> off;
    if (off + a->width > 64) v |= a->w[q + 1] << (64 - off);
    return a->width == 64 ? v : v & ((1ULL << a->width) - 1);
}

#endif
