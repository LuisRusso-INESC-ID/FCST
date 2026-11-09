/* Exhaustive checks of the FCST against a naive suffix tree (SA + LCP). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fcst.h"

static const uint8_t *G;
static size_t        GN;

static int cmp_suf(const void *a, const void *b)
{
    size_t i = *(const size_t *)a, j = *(const size_t *)b;
    while (i < GN && j < GN && G[i] == G[j]) { i++; j++; }
    int ci = i < GN ? G[i] : -1, cj = j < GN ? G[j] : -1;   /* terminator smallest */
    return ci - cj;
}

typedef struct { const uint8_t *T; size_t n; size_t *SA, *ISA, *LCP; } naive_t;

static void naive_build(naive_t *z, const uint8_t *s, size_t m)
{
    z->T = s; z->n = m + 1;
    z->SA = malloc(z->n * sizeof(size_t)); z->ISA = malloc(z->n * sizeof(size_t));
    z->LCP = malloc((z->n + 1) * sizeof(size_t));
    for (size_t i = 0; i < z->n; i++) z->SA[i] = i;
    G = s; GN = m;
    qsort(z->SA, z->n, sizeof(size_t), cmp_suf);
    for (size_t i = 0; i < z->n; i++) z->ISA[z->SA[i]] = i;
    z->LCP[0] = 0;
    for (size_t i = 1; i < z->n; i++) {
        size_t a = z->SA[i - 1], b = z->SA[i], h = 0;
        while (a + h < m && b + h < m && s[a + h] == s[b + h]) h++;
        z->LCP[i] = h;
    }
}

static int nchar(const naive_t *z, size_t p) { return p < z->n - 1 ? z->T[p] : -1; }

static size_t nlcp(const naive_t *z, size_t i, size_t j)
{
    size_t m = (size_t)-1;
    for (size_t k = i + 1; k <= j; k++) if (z->LCP[k] < m) m = z->LCP[k];
    return m;
}

static fcst_node nexpand(const naive_t *z, size_t x, size_t D)
{
    fcst_node v = { x, x };
    while (v.l > 0 && z->LCP[v.l] >= D) v.l--;
    while (v.r + 1 < z->n && z->LCP[v.r + 1] >= D) v.r++;
    return v;
}

static size_t nsdep(const naive_t *z, fcst_node v)
{
    return v.l == v.r ? z->n - z->SA[v.l] : nlcp(z, v.l, v.r);
}

static fcst_node nparent(const naive_t *z, fcst_node v)
{
    if (v.l == 0 && v.r == z->n - 1) return FCST_NONE;
    size_t a = v.l > 0 ? z->LCP[v.l] : 0, b = v.r + 1 < z->n ? z->LCP[v.r + 1] : 0;
    return nexpand(z, v.l, a > b ? a : b);
}

static fcst_node nslink(const naive_t *z, fcst_node v)
{
    if (v.l == 0 && v.r == z->n - 1) return FCST_NONE;
    if (v.l == v.r) {
        size_t p = z->SA[v.l] + 1;
        if (p >= z->n) { fcst_node r = { 0, z->n - 1 }; return r; }
        fcst_node w = { z->ISA[p], z->ISA[p] };
        return w;
    }
    size_t d = nsdep(z, v);
    if (d == 1) { fcst_node r = { 0, z->n - 1 }; return r; }
    return nexpand(z, z->ISA[z->SA[v.l] + 1], d - 1);
}

static long errors = 0;
#define CHECK(cond, ...) do { if (!(cond)) { if (errors++ < 20) { fprintf(stderr, __VA_ARGS__); fputc('\n', stderr); } } } while (0)

static unsigned long long rng = 88172645463325252ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

static void run(const char *name, const uint8_t *s, size_t m, unsigned delta, unsigned rate,
                int fast, unsigned delta_s)
{
    naive_t z;
    naive_build(&z, s, m);
    fcst_t *t = fast ? fcst_build_fast(s, m, delta, delta_s, rate) : fcst_build(s, m, delta, rate);
    size_t  n = z.n;
    long    e0 = errors;

    /* CSA primitives */
    for (size_t i = 0; i < n; i++) {
        CHECK(fcst_locate(t, i) == z.SA[i], "%s locate(%zu)", name, i);
        CHECK(fcst_isa(t, i) == z.ISA[i], "%s isa(%zu)", name, i);
    }

    /* full traversal with FChild / NSibling compared to the naive tree */
    fcst_node *stack = malloc(2 * n * sizeof *stack);
    size_t     top = 0, nodes = 0, internal = 0;
    stack[top++] = fcst_root(t);
    while (top) {
        fcst_node v = stack[--top];
        nodes++;
        size_t d = fcst_sdep(t, v);
        CHECK(d == nsdep(&z, v), "%s sdep [%zu,%zu] %zu vs %zu", name, v.l, v.r, d, nsdep(&z, v));
        fcst_node p = fcst_parent(t, v), np = nparent(&z, v);
        CHECK(fcst_eq(p, np), "%s parent [%zu,%zu]", name, v.l, v.r);
        fcst_node sl = fcst_slink(t, v), nsl = nslink(&z, v);
        CHECK(fcst_eq(sl, nsl), "%s slink [%zu,%zu] got [%zu,%zu] want [%zu,%zu]", name, v.l, v.r, sl.l, sl.r, nsl.l, nsl.r);
        if (d > 0 && !fcst_is_leaf(v)) {
            size_t k = 1 + xr() % (d + 2);
            fcst_node w = fcst_slink_i(t, v, k), nw = v;
            for (size_t q = 0; q < k && !(nw.l == 0 && nw.r == n - 1); q++) nw = nslink(&z, nw);
            CHECK(fcst_eq(w, nw), "%s slink_i [%zu,%zu] k=%zu", name, v.l, v.r, k);
            size_t i = 1 + xr() % d;
            int c = fcst_letter(t, v, i);
            CHECK(c == nchar(&z, z.SA[v.l] + i - 1), "%s letter", name);
        }
        if (fcst_is_leaf(v)) continue;
        internal++;
        /* naive children */
        size_t    nc = 0, x = v.l;
        fcst_node ch = fcst_fchild(t, v);
        while (x <= v.r) {
            int    c = nchar(&z, z.SA[x] + d);
            size_t y = x;
            while (y + 1 <= v.r && nchar(&z, z.SA[y + 1] + d) == c) y++;
            fcst_node want = { x, y };
            CHECK(!fcst_is_none(ch) && fcst_eq(ch, want), "%s child #%zu of [%zu,%zu]", name, nc, v.l, v.r);
            if (c >= 0) CHECK(fcst_eq(fcst_child(t, v, (uint8_t)c), want), "%s child(c)", name);
            if (fcst_is_none(ch)) break;
            stack[top++] = ch;
            ch = fcst_nsibling(t, ch);
            nc++;
            x = y + 1;
        }
        CHECK(fcst_is_none(ch), "%s extra sibling under [%zu,%zu]", name, v.l, v.r);
    }
    free(stack);
    CHECK(nodes == n + internal, "%s node count", name);

    /* random LCAs of nodes */
    for (int it = 0; it < 2000; it++) {
        size_t i = xr() % n, j = xr() % n;
        if (i > j) { size_t q = i; i = j; j = q; }
        fcst_node a = { i, i }, b = { j, j };
        fcst_node got = fcst_lca(t, a, b);
        fcst_node want = i == j ? a : nexpand(&z, i, nlcp(&z, i, j));
        CHECK(fcst_eq(got, want), "%s lca(%zu,%zu)", name, i, j);
        fcst_node pa = fcst_parent(t, a), pb = fcst_parent(t, b);
        if (fcst_is_none(pa) || fcst_is_none(pb)) continue;
        fcst_node g2 = fcst_lca(t, pa, pb);
        fcst_node w2 = nexpand(&z, pa.l < pb.l ? pa.l : pb.l,
                               nlcp(&z, pa.l < pb.l ? pa.l : pb.l, pa.r > pb.r ? pa.r : pb.r));
        if (fcst_ancestor(pa, pb)) w2 = pa; else if (fcst_ancestor(pb, pa)) w2 = pb;
        CHECK(fcst_eq(g2, w2), "%s lca(parents)", name);
    }

    /* pattern search and Weiner links */
    for (int it = 0; it < 300 && m > 0; it++) {
        size_t p = xr() % m, len = 1 + xr() % 8;
        if (p + len > m) len = m - p;
        fcst_node v = fcst_search(t, s + p, len);
        size_t occ = 0;
        for (size_t q = 0; q + len <= m; q++) occ += memcmp(s + q, s + p, len) == 0;
        CHECK(!fcst_is_none(v) && fcst_count(v) == occ, "%s search", name);
        if (p > 0) {
            fcst_node w = fcst_search(t, s + p, len);
            fcst_node wl = fcst_wlink(t, w, s[p - 1]), want = fcst_search(t, s + p - 1, len + 1);
            CHECK(fcst_eq(wl, want), "%s wlink", name);
        }
    }

    printf("%-5s %-12s n=%-6zu delta=%-3u", fast ? "fast" : "orig", name, n, delta);
    if (fast) printf(" delta_s=%-3u", delta_s); else printf("            ");
    printf(" sampled=%-6zu internal=%-6zu %s\n", fcst_sampled(t), internal, errors == e0 ? "OK" : "FAIL");
    fcst_free(t);
    free(z.SA); free(z.ISA); free(z.LCP);
}

int main(void)
{
    size_t   M = 3000;
    uint8_t *s = malloc(M + 1);
    unsigned deltas[] = { 2, 3, 4, 6, 8, 16, 33, 64 };
    unsigned dsdeps[] = { 8, 2, 16, 32, 64, 4, 64, 128 };   /* paired with deltas in fast mode */

    for (size_t ci = 0; ci < 2 * sizeof deltas / sizeof *deltas; ci++) {
        size_t   di = ci % (sizeof deltas / sizeof *deltas);
        int      fast = ci >= sizeof deltas / sizeof *deltas;
        unsigned d = deltas[di], ds = dsdeps[di];
        size_t   m;
        for (m = 0; m < 2000; m++) s[m] = "ab"[xr() % 2];
        run("random-2", s, m, d, 4 + d % 7, fast, ds);
        for (m = 0; m < 2000; m++) s[m] = "ACGT"[xr() % 4];
        run("dna", s, m, d, 16, fast, ds);
        for (m = 0; m < 1500; m++) s[m] = (uint8_t)('a' + xr() % 26);
        run("random-26", s, m, d, 8, fast, ds);
        for (m = 0; m < 1200; m++) s[m] = 'a';
        run("aaaa", s, m, d, 32, fast, ds);
        for (m = 0; m < 1200; m++) s[m] = "abaab"[m % 5];
        run("periodic", s, m, d, 7, fast, ds);
        {   /* Fibonacci word */
            char *a = calloc(4000, 1), *b = calloc(4000, 1), *c = calloc(4000, 1);
            strcpy(a, "a"); strcpy(b, "ab");
            while (strlen(b) < 1500) { strcpy(c, b); strcat(c, a); strcpy(a, b); strcpy(b, c); }
            m = strlen(b); memcpy(s, b, m);
            free(a); free(b); free(c);
        }
        run("fibonacci", s, m, d, 5, fast, ds);
        m = 1; s[0] = 'x';
        run("single", s, m, d, 3, fast, ds);
        m = 0;
        run("empty", s, m, d, 3, fast, ds);
        for (m = 0; m < 2500; m++) s[m] = (uint8_t)(1 + xr() % 255);
        run("bytes-255", s, m, d, 9, fast, ds);
    }
    free(s);
    printf(errors ? "\n%ld ERRORS\n" : "\nALL TESTS PASSED\n", errors);
    return errors != 0;
}
