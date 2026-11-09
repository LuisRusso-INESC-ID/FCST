#include "bits.h"

#include <stdlib.h>
#include <string.h>

#define SB_WORDS 8          /* 512-bit superblocks          */
#define SEL      4096       /* select sampling rate         */

static inline unsigned popc(uint64_t x) { return (unsigned)__builtin_popcountll(x); }

static inline unsigned select_in_word(uint64_t x, size_t r) /* r >= 1 */
{
    while (--r) x &= x - 1;
    return (unsigned)__builtin_ctzll(x);
}

void bv_init(bv_t *b, size_t n)
{
    memset(b, 0, sizeof *b);
    b->n  = n;
    b->nw = (n + 63) / 64;
    b->w  = calloc(b->nw + 1, sizeof(uint64_t));
}

void bv_build(bv_t *b)
{
    size_t nsb = (b->nw + SB_WORDS - 1) / SB_WORDS;
    if (nsb == 0) nsb = 1;
    b->nsb = nsb;
    b->sb  = malloc((nsb + 1) * sizeof(uint64_t));
    uint64_t c = 0;
    for (size_t j = 0; j < nsb; j++) {
        b->sb[j] = c;
        for (size_t q = j * SB_WORDS; q < (j + 1) * SB_WORDS && q < b->nw; q++) c += popc(b->w[q]);
    }
    b->sb[nsb] = c;
    b->ones    = c;

    size_t zeros = b->n - c;
    b->ns1 = c / SEL + 1;
    b->ns0 = zeros / SEL + 1;
    b->s1  = malloc(b->ns1 * sizeof(uint32_t));
    b->s0  = malloc(b->ns0 * sizeof(uint32_t));
    size_t j = 0;
    for (size_t m = 0; m < b->ns1; m++) {
        uint64_t t = (uint64_t)m * SEL + 1;
        while (j + 1 < nsb && b->sb[j + 1] < t) j++;
        b->s1[m] = (uint32_t)j;
    }
    j = 0;
    for (size_t m = 0; m < b->ns0; m++) {
        uint64_t t = (uint64_t)m * SEL + 1;
        while (j + 1 < nsb && (j + 1) * 512 - b->sb[j + 1] < t) j++;
        b->s0[m] = (uint32_t)j;
    }
}

void bv_free(bv_t *b)
{
    free(b->w); free(b->sb); free(b->s1); free(b->s0);
    memset(b, 0, sizeof *b);
}

size_t bv_rank1(const bv_t *b, size_t i)
{
    size_t j = i / 512, c = b->sb[j], end = i >> 6;
    for (size_t q = j * SB_WORDS; q < end; q++) c += popc(b->w[q]);
    if (i & 63) c += popc(b->w[end] & ((1ULL << (i & 63)) - 1));
    return c;
}

size_t bv_select1(const bv_t *b, size_t k)
{
    size_t j = b->s1[(k - 1) / SEL];
    while (j + 1 < b->nsb && b->sb[j + 1] < k) j++;
    size_t r = k - b->sb[j], q = j * SB_WORDS;
    for (;; q++) {
        unsigned pc = popc(b->w[q]);
        if (r <= pc) break;
        r -= pc;
    }
    return q * 64 + select_in_word(b->w[q], r);
}

size_t bv_select0(const bv_t *b, size_t k)
{
    size_t j = b->s0[(k - 1) / SEL];
    while (j + 1 < b->nsb && (j + 1) * 512 - b->sb[j + 1] < k) j++;
    size_t r = k - (j * 512 - b->sb[j]), q = j * SB_WORDS;
    for (;; q++) {
        unsigned pc = popc(~b->w[q]);
        if (r <= pc) break;
        r -= pc;
    }
    return q * 64 + select_in_word(~b->w[q], r);
}

size_t bv_bytes(const bv_t *b)
{
    return (b->nw + 1) * 8 + (b->nsb + 1) * 8 + (b->ns1 + b->ns0) * 4;
}

unsigned bits_needed(uint64_t x)
{
    unsigned w = 1;
    while (w < 64 && (x >> w)) w++;
    return w;
}

void pa_init(pa_t *a, size_t n, unsigned width)
{
    a->n     = n;
    a->width = width ? width : 1;
    a->w     = calloc((n * a->width + 63) / 64 + 1, sizeof(uint64_t));
}

void pa_free(pa_t *a) { free(a->w); a->w = NULL; }

void pa_set(pa_t *a, size_t i, uint64_t v)
{
    size_t   pos = i * a->width, q = pos >> 6;
    unsigned off = pos & 63;
    uint64_t mask = a->width == 64 ? ~0ULL : (1ULL << a->width) - 1;
    v &= mask;
    a->w[q] = (a->w[q] & ~(mask << off)) | (v << off);
    if (off + a->width > 64) {
        unsigned hi = off + a->width - 64;
        uint64_t m2 = (1ULL << hi) - 1;
        a->w[q + 1] = (a->w[q + 1] & ~m2) | (v >> (64 - off));
    }
}

size_t pa_bytes(const pa_t *a) { return ((a->n * a->width + 63) / 64 + 1) * 8; }
