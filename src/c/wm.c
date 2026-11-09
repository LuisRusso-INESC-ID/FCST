#include "wm.h"

#include <stdlib.h>
#include <string.h>

void wm_build(wm_t *w, const uint8_t *s, size_t n, unsigned sigma)
{
    w->n  = n;
    w->L  = bits_needed(sigma > 1 ? sigma - 1 : 1);
    w->lv = calloc(w->L, sizeof(bv_t));
    w->z  = calloc(w->L, sizeof(size_t));
    uint8_t *cur = malloc(n ? n : 1), *nxt = malloc(n ? n : 1);
    memcpy(cur, s, n);
    for (unsigned l = 0; l < w->L; l++) {
        unsigned sh = w->L - 1 - l;
        bv_init(&w->lv[l], n);
        size_t zc = 0;
        for (size_t i = 0; i < n; i++) {
            if ((cur[i] >> sh) & 1) bv_set(&w->lv[l], i);
            else zc++;
        }
        w->z[l] = zc;
        size_t a = 0, b = zc;
        for (size_t i = 0; i < n; i++) {
            if ((cur[i] >> sh) & 1) nxt[b++] = cur[i];
            else nxt[a++] = cur[i];
        }
        uint8_t *t = cur; cur = nxt; nxt = t;
        bv_build(&w->lv[l]);
    }
    free(cur); free(nxt);
}

void wm_free(wm_t *w)
{
    for (unsigned l = 0; l < w->L; l++) bv_free(&w->lv[l]);
    free(w->lv); free(w->z);
}

unsigned wm_access(const wm_t *w, size_t i)
{
    unsigned c = 0;
    for (unsigned l = 0; l < w->L; l++) {
        const bv_t *b = &w->lv[l];
        int bit = bv_get(b, i);
        c = (c << 1) | (unsigned)bit;
        i = bit ? w->z[l] + bv_rank1(b, i) : bv_rank0(b, i);
    }
    return c;
}

size_t wm_rank(const wm_t *w, unsigned c, size_t i)
{
    size_t s = 0;
    for (unsigned l = 0; l < w->L; l++) {
        const bv_t *b = &w->lv[l];
        if ((c >> (w->L - 1 - l)) & 1) {
            i = w->z[l] + bv_rank1(b, i);
            s = w->z[l] + bv_rank1(b, s);
        } else {
            i = bv_rank0(b, i);
            s = bv_rank0(b, s);
        }
    }
    return i - s;
}

unsigned wm_access_rank(const wm_t *w, size_t i, size_t *rank)
{
    unsigned c = 0;
    size_t   s = 0;
    for (unsigned l = 0; l < w->L; l++) {
        const bv_t *b = &w->lv[l];
        int bit = bv_get(b, i);
        c = (c << 1) | (unsigned)bit;
        if (bit) {
            i = w->z[l] + bv_rank1(b, i);
            s = w->z[l] + bv_rank1(b, s);
        } else {
            i = bv_rank0(b, i);
            s = bv_rank0(b, s);
        }
    }
    *rank = i - s;
    return c;
}

size_t wm_select(const wm_t *w, unsigned c, size_t k)
{
    size_t s = 0;
    for (unsigned l = 0; l < w->L; l++) {
        const bv_t *b = &w->lv[l];
        s = ((c >> (w->L - 1 - l)) & 1) ? w->z[l] + bv_rank1(b, s) : bv_rank0(b, s);
    }
    size_t p = s + k - 1;
    for (unsigned l = w->L; l-- > 0;) {
        const bv_t *b = &w->lv[l];
        if ((c >> (w->L - 1 - l)) & 1) p = bv_select1(b, p - w->z[l] + 1);
        else                           p = bv_select0(b, p + 1);
    }
    return p;
}

size_t wm_bytes(const wm_t *w)
{
    size_t t = sizeof *w + w->L * sizeof(size_t);
    for (unsigned l = 0; l < w->L; l++) t += bv_bytes(&w->lv[l]);
    return t;
}
