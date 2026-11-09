/* usage: bench FILE [spec ...]  spec = delta (original FCST) or dpar:dsdep (fast FCST)
 * sa_rate via env FCST_SA_RATE, default 32 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "fcst.h"

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static unsigned long long rng = 0x9E3779B97F4A7C15ULL;
static uint64_t xr(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng; }

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: %s FILE [delta ...]\n", argv[0]); return 1; }
    FILE *f = fopen(argv[1], "rb");
    if (!f) { perror(argv[1]); return 1; }
    fseek(f, 0, SEEK_END);
    size_t m = (size_t)ftell(f);
    rewind(f);
    uint8_t *s = malloc(m + 1);
    if (fread(s, 1, m, f) != m) { perror("read"); return 1; }
    fclose(f);
    for (size_t i = 0; i < m; i++) if (!s[i]) s[i] = 1;   /* no zero bytes */

    unsigned rate = getenv("FCST_SA_RATE") ? (unsigned)atoi(getenv("FCST_SA_RATE")) : 32;
    int      nd = argc > 2 ? argc - 2 : 1;
    const int Q = 20000;

    for (int di = 0; di < nd; di++) {
        const char *spec = argc > 2 ? argv[2 + di] : "8";
        unsigned delta = (unsigned)atoi(spec), dsd = 0;
        const char *colon = strchr(spec, ':');
        if (colon) dsd = (unsigned)atoi(colon + 1);
        double   t0 = now();
        fcst_t  *t = colon ? fcst_build_fast(s, m, delta, dsd, rate) : fcst_build(s, m, delta, rate);
        double   tb = now() - t0;
        size_t   n = fcst_n(t);
        printf("\n== %s  %s %s  sa_rate=%u  build %.2fs\n", argv[1], colon ? "fast" : "orig", spec, rate, tb);
        fcst_space(t, stdout);

        /* sample internal nodes as parents of random leaves */
        fcst_node *V = malloc(Q * sizeof *V);
        rng = 12345;
        t0 = now();
        for (int q = 0; q < Q; q++) {
            size_t i = xr() % n;
            fcst_node leaf = { i, i };
            V[q] = fcst_parent(t, leaf);
        }
        double tp = (now() - t0) / Q;
        size_t chk = 0;
        t0 = now();
        for (int q = 0; q < Q; q++) chk += fcst_sdep(t, V[q]);
        double ts = (now() - t0) / Q;
        t0 = now();
        for (int q = 0; q < Q; q++) chk += fcst_slink(t, V[q]).l;
        double tl = (now() - t0) / Q;
        t0 = now();
        for (int q = 0; q + 1 < Q; q++) chk += fcst_lca(t, V[q], V[q + 1]).l;
        double ta = (now() - t0) / Q;
        t0 = now();
        for (int q = 0; q < Q / 10; q++) chk += fcst_fchild(t, V[q]).l;
        double tc = (now() - t0) / (Q / 10);
        t0 = now();
        for (int q = 0; q < Q; q++) chk += fcst_locate(t, xr() % n);
        double tloc = (now() - t0) / Q;
        t0 = now();
        for (int q = 0; q < Q; q++) chk += fcst_psi(t, xr() % n);
        double tpsi = (now() - t0) / Q;
        printf("  us/op: Psi %.2f  Locate %.2f  Parent %.2f  SDep %.2f  SLink %.2f  LCA %.2f  FChild %.2f  [%zu]\n",
               tpsi * 1e6, tloc * 1e6, tp * 1e6, ts * 1e6, tl * 1e6, ta * 1e6, tc * 1e6, chk % 10);
        free(V);
        fcst_free(t);
    }
    free(s);
    return 0;
}
