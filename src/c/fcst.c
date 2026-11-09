#include "fcst.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

#include "bits.h"
#include "wm.h"

struct fcst {
    size_t   n;              /* |T$|                                   */
    unsigned sigma;          /* alphabet size including terminator      */
    uint8_t  sym[256];       /* code -> byte                            */
    uint8_t  map[256];       /* byte -> code (0 = absent)               */
    size_t   C[257];         /* C[c] = # symbols < c                    */
    wm_t     bwt;

    unsigned sa_rate;
    bv_t     sa_mark;        /* SA positions whose A[i] = 0 mod rate    */
    pa_t     sa_samp;        /* A[i] / rate for marked i                */
    pa_t     isa_samp;       /* A^{-1}[k*rate]                          */

    unsigned h;              /* delta / 2                               */
    size_t   ns;             /* sampled nodes                           */
    bv_t     B;              /* leaves (0) interleaved with parens (1)  */
    bv_t     P;              /* parens: 1 = '(' 0 = ')'                 */
    size_t   np;             /* 2*ns                                    */
    size_t   nb, segsz;      /* blocks of 64 excess positions           */
    int32_t *seg;            /* min-segment tree over block minima      */
    pa_t     sd;             /* original mode: SDep / h of every sampled node */

    /* fast mode: LCA by smallest interval; SDep stored only on a sparser set */
    int       fast;
    unsigned  hs;            /* delta_sdep / 2                            */
    size_t    nsd;           /* nodes with a stored SDep                  */
    bv_t      M;             /* over sampled nodes (preorder): SDep stored */
    uint64_t *ds;            /* Elias-delta codes of SDep/hs + 1          */
    size_t    dsbits;
    pa_t      dptr;          /* bit offset of every DBLK-th code          */
};

#define DBLK 16

/* ------------------------------------------------------------------ */
/* FM-index primitives                                                 */

static inline unsigned fcode(const fcst_t *t, size_t i)   /* F[i] */
{
    unsigned lo = 0, hi = t->sigma - 1;
    while (lo < hi) {
        unsigned mid = (lo + hi + 1) / 2;
        if (t->C[mid] <= i) lo = mid; else hi = mid - 1;
    }
    return lo;
}

size_t fcst_psi(const fcst_t *t, size_t i)
{
    unsigned c = fcode(t, i);
    return wm_select(&t->bwt, c, i - t->C[c] + 1);
}

size_t fcst_lf(const fcst_t *t, size_t i)
{
    size_t   r;
    unsigned c = wm_access_rank(&t->bwt, i, &r);
    return t->C[c] + r;
}

static inline int lf_interval(const fcst_t *t, unsigned c, size_t *l, size_t *r)
{
    size_t nl = t->C[c] + wm_rank(&t->bwt, c, *l);
    size_t nr = t->C[c] + wm_rank(&t->bwt, c, *r + 1);
    if (nl >= nr) return 0;
    *l = nl; *r = nr - 1;
    return 1;
}

size_t fcst_locate(const fcst_t *t, size_t i)
{
    size_t steps = 0;
    while (!bv_get(&t->sa_mark, i)) { i = fcst_lf(t, i); steps++; }
    return pa_get(&t->sa_samp, bv_rank1(&t->sa_mark, i)) * t->sa_rate + steps;
}

size_t fcst_isa(const fcst_t *t, size_t p)
{
    if (p >= t->n - 1) return 0;                 /* the terminator suffix */
    size_t tp = (p + t->sa_rate - 1) / t->sa_rate * t->sa_rate, i;
    if (tp >= t->n - 1) { tp = t->n - 1; i = 0; }
    else i = pa_get(&t->isa_samp, tp / t->sa_rate);
    for (; tp > p; tp--) i = fcst_lf(t, i);
    return i;
}

/* code of T[A[x] + d] */
static unsigned code_at(const fcst_t *t, size_t x, size_t d)
{
    if (d <= t->sa_rate) {
        while (d--) x = fcst_psi(t, x);
        return fcode(t, x);
    }
    return fcode(t, fcst_isa(t, fcst_locate(t, x) + d));
}

/* ------------------------------------------------------------------ */
/* Sampled tree: parentheses, excess, LCSA                             */

static inline long excess(const fcst_t *t, size_t k)
{
    return 2 * (long)bv_rank1(&t->P, k) - (long)k;
}

static long scan_min(const fcst_t *t, size_t k0, size_t k1)  /* min E(k), k in [k0,k1] */
{
    long e = excess(t, k0), m = e;
    for (size_t k = k0; k < k1; k++) {
        e += bv_get(&t->P, k) ? 1 : -1;
        if (e < m) m = e;
    }
    return m;
}

static long seg_min(const fcst_t *t, size_t L, size_t R)   /* blocks [L,R] */
{
    long m = LONG_MAX;
    L += t->segsz; R += t->segsz + 1;
    while (L < R) {
        if (L & 1) { if (t->seg[L] < m) m = t->seg[L]; L++; }
        if (R & 1) { --R; if (t->seg[R] < m) m = t->seg[R]; }
        L >>= 1; R >>= 1;
    }
    return m;
}

/* rightmost block index < R whose minimum is <= th, or -1 */
static long seg_rightmost(const fcst_t *t, size_t node, size_t nl, size_t nr, size_t R, long th)
{
    if (nl >= R || t->seg[node] > th) return -1;
    if (nr - nl == 1) return (long)nl;
    size_t mid = (nl + nr) / 2;
    long   r   = seg_rightmost(t, 2 * node + 1, mid, nr, R, th);
    return r >= 0 ? r : seg_rightmost(t, 2 * node, nl, mid, R, th);
}

static long rmq_min(const fcst_t *t, size_t a, size_t b)
{
    size_t ba = a / 64, bb = b / 64;
    if (ba == bb) return scan_min(t, a, b);
    long m = scan_min(t, a, ba * 64 + 63), x;
    if (ba + 1 <= bb - 1 && (x = seg_min(t, ba + 1, bb - 1)) < m) m = x;
    if ((x = scan_min(t, bb * 64, b)) < m) m = x;
    return m;
}

/* largest p < a with E(p) <= th */
static size_t bwd_search(const fcst_t *t, size_t a, long th)
{
    size_t p = a, start = (a / 64) * 64;
    long   e = excess(t, a);
    while (p > start) {
        p--;
        e -= bv_get(&t->P, p) ? 1 : -1;
        if (e <= th) return p;
    }
    long b = seg_rightmost(t, 1, 0, t->segsz, a / 64, th);
    size_t k0 = (size_t)b * 64, k1 = k0 + 63, best = k0;
    if (k1 > t->np) k1 = t->np;
    e = excess(t, k0);
    for (size_t k = k0; k < k1; k++) {
        e += bv_get(&t->P, k) ? 1 : -1;
        if (e <= th) best = k + 1;
    }
    return best;
}

/* leftmost block index >= L whose minimum is <= th, or -1 */
static long seg_leftmost(const fcst_t *t, size_t node, size_t nl, size_t nr, size_t L, long th)
{
    if (nr <= L || t->seg[node] > th) return -1;
    if (nr - nl == 1) return (long)nl;
    size_t mid = (nl + nr) / 2;
    long   r   = seg_leftmost(t, 2 * node, nl, mid, L, th);
    return r >= 0 ? r : seg_leftmost(t, 2 * node + 1, mid, nr, L, th);
}

/* smallest k > p with E(k) <= th */
static size_t fwd_search(const fcst_t *t, size_t p, long th)
{
    size_t k = p, end = (p / 64) * 64 + 63;
    long   e = excess(t, p);
    if (end > t->np) end = t->np;
    while (k < end) {
        e += bv_get(&t->P, k) ? 1 : -1;
        k++;
        if (e <= th) return k;
    }
    long b = seg_leftmost(t, 1, 0, t->segsz, p / 64 + 1, th);
    k = (size_t)b * 64;
    e = excess(t, k);
    while (e > th) { e += bv_get(&t->P, k) ? 1 : -1; k++; }
    return k;
}

/* SA interval of the sampled node whose '(' is P[p] */
static void sampled_interval(const fcst_t *t, size_t p, size_t *l, size_t *r)
{
    size_t q = fwd_search(t, p, excess(t, p)) - 1;      /* matching ')' */
    *l = bv_select1(&t->B, p + 1) - p;
    *r = bv_select1(&t->B, q + 1) - q - 1;
}

/* lowest common sampled ancestor of leaves i < j: position of its '(' in P */
static size_t lcsa(const fcst_t *t, size_t i, size_t j)
{
    size_t a = bv_select0(&t->B, i + 1) - i;
    size_t b = bv_select0(&t->B, j + 1) - j;
    long   m = rmq_min(t, a, b);
    return bwd_search(t, a, m - 1);
}

/* ------------------------------------------------------------------ */
/* Elias-delta codes (MSB first)                                       */

static void put_bits(uint64_t *w, size_t *pos, uint64_t v, unsigned nb)
{
    for (unsigned q = nb; q-- > 0; (*pos)++)
        if ((v >> q) & 1) w[*pos >> 6] |= 1ULL << (*pos & 63);
}

static unsigned delta_len(uint64_t x)               /* x >= 1 */
{
    unsigned N = 63 - (unsigned)__builtin_clzll(x), L = N + 1;
    unsigned LL = 31 - (unsigned)__builtin_clz(L);
    return 2 * LL + 1 + N;
}

static void put_delta(uint64_t *w, size_t *pos, uint64_t x)
{
    unsigned N = 63 - (unsigned)__builtin_clzll(x), L = N + 1;
    unsigned LL = 31 - (unsigned)__builtin_clz(L);
    *pos += LL;                                       /* LL zeros */
    put_bits(w, pos, L, LL + 1);
    put_bits(w, pos, x, N);                           /* low N bits */
}

static inline int get_bit(const uint64_t *w, size_t pos) { return (int)((w[pos >> 6] >> (pos & 63)) & 1); }

static uint64_t get_delta(const uint64_t *w, size_t *pos)
{
    unsigned LL = 0;
    while (!get_bit(w, *pos)) { LL++; (*pos)++; }
    uint64_t L = 0;
    for (unsigned q = 0; q <= LL; q++) L = (L << 1) | (uint64_t)get_bit(w, (*pos)++);
    uint64_t x = 1;
    for (uint64_t q = 1; q < L; q++) x = (x << 1) | (uint64_t)get_bit(w, (*pos)++);
    return x;
}

/* stored SDep of the sampled node with preorder id, or -1 */
static long stored_sdep(const fcst_t *t, size_t id)
{
    if (!t->fast) return (long)pa_get(&t->sd, id) * t->h;
    if (!bv_get(&t->M, id)) return -1;
    size_t r = bv_rank1(&t->M, id), pos = pa_get(&t->dptr, r / DBLK);
    for (size_t q = r % DBLK; q > 0; q--) get_delta(t->ds, &pos);
    return (long)(get_delta(t->ds, &pos) - 1) * t->hs;
}

/* ------------------------------------------------------------------ */
/* Fast mode                                                           */

/* LCA of leaves i < j without string depths: for each k the LCSA mapped
 * back with LF is an ancestor of the LCA, and the LCA is the smallest
 * such interval.  B_k = min(LCSA_k, LF(letter_k, B_{k+1})).               */
static fcst_node lca_fast(const fcst_t *t, size_t i, size_t j)
{
    unsigned lim = 2 * t->h;
    uint8_t  lbuf[256], *let = lim <= sizeof lbuf ? lbuf : malloc(lim);
    size_t   pbuf[256], *par = lim <= 256 ? pbuf : malloc(lim * sizeof *par);
    unsigned kk = 0;
    for (unsigned k = 0; k < lim; k++) {
        unsigned ci = fcode(t, i), cj = fcode(t, j);
        kk = k + 1;
        if (ci != cj) { par[k] = 0; break; }          /* root */
        par[k] = lcsa(t, i, j);
        let[k] = (uint8_t)ci;
        i = fcst_psi(t, i);
        j = fcst_psi(t, j);
    }
    size_t l, r, l2, r2;
    sampled_interval(t, par[kk - 1], &l, &r);
    for (unsigned k = kk - 1; k-- > 0;) {
        lf_interval(t, let[k], &l, &r);
        if (par[k] == 0) continue;                     /* root never smaller */
        sampled_interval(t, par[k], &l2, &r2);
        if (r2 - l2 < r - l) { l = l2; r = r2; }
    }
    if (let != lbuf) free(let);
    if (par != pbuf) free(par);
    fcst_node v = { l, r };
    return v;
}

/* SDep(LCA(i,j)), i < j: max over k < 2hs of k + stored SDep(LCSA_k) */
static size_t sdep_fast(const fcst_t *t, size_t i, size_t j)
{
    unsigned lim = 2 * t->hs;
    size_t   best = 0;
    for (unsigned k = 0; k < lim; k++) {
        if (fcode(t, i) != fcode(t, j)) { if (k > best) best = k; break; }
        size_t p = lcsa(t, i, j);
        long   s = stored_sdep(t, bv_rank1(&t->P, p));
        if (s >= 0 && k + (size_t)s > best) best = k + (size_t)s;
        i = fcst_psi(t, i);
        j = fcst_psi(t, j);
    }
    return best;
}

/* ------------------------------------------------------------------ */
/* LCA of two leaves i < j: returns SDep, optionally the node          */

static size_t lca_leaves(const fcst_t *t, size_t i, size_t j, fcst_node *out)
{
    unsigned lim = 2 * t->h;
    uint8_t  buf[512], *let = lim <= sizeof buf ? buf : malloc(lim);
    size_t   best = 0, bestk = 0, bestp = 0;   /* bestp: '(' of best LCSA */

    for (unsigned k = 0; k < lim; k++) {
        unsigned ci = fcode(t, i), cj = fcode(t, j);
        if (ci != cj) {                      /* LCSA is the root here */
            if (k >= best) { best = k; bestk = k; bestp = 0; }
            break;
        }
        size_t p    = lcsa(t, i, j);
        size_t cand = k + (size_t)pa_get(&t->sd, bv_rank1(&t->P, p)) * t->h;
        if (cand > best || k == 0) { best = cand; bestk = k; bestp = p; }
        let[k] = (uint8_t)ci;
        i = fcst_psi(t, i);
        j = fcst_psi(t, j);
    }
    if (out) {
        size_t l, r;
        sampled_interval(t, bestp, &l, &r);
        for (size_t q = bestk; q-- > 0;) lf_interval(t, let[q], &l, &r);
        out->l = l; out->r = r;
    }
    if (let != buf) free(let);
    return best;
}

static fcst_node lca_node(const fcst_t *t, size_t i, size_t j)
{
    fcst_node v;
    if (t->fast) return lca_fast(t, i, j);
    lca_leaves(t, i, j, &v);
    return v;
}

/* ------------------------------------------------------------------ */
/* Tree operations                                                     */

size_t    fcst_n(const fcst_t *t)       { return t->n; }
size_t    fcst_sampled(const fcst_t *t) { return t->ns; }
fcst_node fcst_root(const fcst_t *t)    { fcst_node v = { 0, t->n - 1 }; return v; }

static inline int is_root(const fcst_t *t, fcst_node v) { return v.l == 0 && v.r == t->n - 1; }

size_t fcst_sdep(const fcst_t *t, fcst_node v)
{
    if (fcst_is_leaf(v)) return t->n - fcst_locate(t, v.l);
    return t->fast ? sdep_fast(t, v.l, v.r) : lca_leaves(t, v.l, v.r, NULL);
}

fcst_node fcst_lca(const fcst_t *t, fcst_node v, fcst_node w)
{
    if (fcst_ancestor(v, w)) return v;
    if (fcst_ancestor(w, v)) return w;
    return lca_node(t, v.l < w.l ? v.l : w.l, v.r > w.r ? v.r : w.r);
}

fcst_node fcst_parent(const fcst_t *t, fcst_node v)
{
    if (is_root(t, v)) return FCST_NONE;
    fcst_node best = FCST_NONE;
    if (v.l > 0) best = lca_node(t, v.l - 1, v.l);
    if (v.r < t->n - 1) {
        fcst_node y = lca_node(t, v.r, v.r + 1);
        if (fcst_is_none(best) || fcst_count(y) < fcst_count(best)) best = y;
    }
    return best;
}

fcst_node fcst_slink(const fcst_t *t, fcst_node v)
{
    if (is_root(t, v)) return FCST_NONE;
    if (fcst_is_leaf(v)) {
        if (v.l == 0) return fcst_root(t);          /* suffix "$" */
        size_t x = fcst_psi(t, v.l);
        fcst_node w = { x, x };
        return w;
    }
    return lca_node(t, fcst_psi(t, v.l), fcst_psi(t, v.r));
}

fcst_node fcst_slink_i(const fcst_t *t, fcst_node v, size_t k)
{
    if (k == 0) return v;
    if (is_root(t, v)) return FCST_NONE;
    if (fcst_is_leaf(v)) {
        size_t x = v.l;
        for (size_t q = 0; q < k; q++) {
            if (x == 0) return fcst_root(t);
            x = fcst_psi(t, x);
        }
        fcst_node w = { x, x };
        return w;
    }
    size_t x = v.l, y = v.r;
    for (size_t q = 0; q < k; q++) {
        if (fcode(t, x) != fcode(t, y)) return fcst_root(t);   /* k >= SDep(v) */
        x = fcst_psi(t, x);
        y = fcst_psi(t, y);
    }
    return lca_node(t, x, y);
}

static fcst_node child_code(const fcst_t *t, fcst_node v, size_t d, unsigned c)
{
    size_t lo = v.l, hi = v.r + 1;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (code_at(t, mid, d) < c) lo = mid + 1; else hi = mid;
    }
    size_t a = lo;
    hi = v.r + 1;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (code_at(t, mid, d) <= c) lo = mid + 1; else hi = mid;
    }
    if (a >= lo) return FCST_NONE;
    fcst_node w = { a, lo - 1 };
    return w;
}

fcst_node fcst_child(const fcst_t *t, fcst_node v, uint8_t c)
{
    if (fcst_is_leaf(v) || c == 0 || !t->map[c]) return FCST_NONE;
    return child_code(t, v, fcst_sdep(t, v), t->map[c]);
}

fcst_node fcst_fchild(const fcst_t *t, fcst_node v)
{
    if (fcst_is_leaf(v)) return FCST_NONE;
    size_t d = fcst_sdep(t, v);
    return child_code(t, v, d, code_at(t, v.l, d));
}

fcst_node fcst_nsibling(const fcst_t *t, fcst_node v)
{
    if (is_root(t, v)) return FCST_NONE;
    fcst_node p = fcst_parent(t, v);
    if (v.r == p.r) return FCST_NONE;
    size_t d = fcst_sdep(t, p);
    return child_code(t, p, d, code_at(t, v.r + 1, d));
}

int fcst_letter(const fcst_t *t, fcst_node v, size_t i)
{
    unsigned c = code_at(t, v.l, i - 1);
    return c == 0 ? -1 : t->sym[c];
}

fcst_node fcst_wlink(const fcst_t *t, fcst_node v, uint8_t c)
{
    if (c == 0 || !t->map[c]) return FCST_NONE;
    size_t l = v.l, r = v.r;
    if (!lf_interval(t, t->map[c], &l, &r)) return FCST_NONE;
    fcst_node w = { l, r };
    return w;
}

fcst_node fcst_search(const fcst_t *t, const uint8_t *p, size_t m)
{
    size_t l = 0, r = t->n - 1;
    for (size_t q = m; q-- > 0;) {
        if (p[q] == 0 || !t->map[p[q]]) return FCST_NONE;
        if (!lf_interval(t, t->map[p[q]], &l, &r)) return FCST_NONE;
    }
    fcst_node w = { l, r };
    return w;
}

/* ------------------------------------------------------------------ */
/* Construction (uses O(n log n) bits temporarily)                     */

static uint32_t *build_sa(const uint8_t *T, size_t n, unsigned sigma)
{
    uint32_t *sa  = malloc(n * sizeof *sa), *rk = malloc(n * sizeof *rk);
    uint32_t *tmp = malloc(n * sizeof *tmp);
    size_t    csz = (n > sigma ? n : sigma) + 2;
    uint32_t *cnt = malloc(csz * sizeof *cnt);
    for (size_t i = 0; i < n; i++) rk[i] = T[i];
    size_t maxr = sigma - 1;
    if (n == 1) { sa[0] = 0; goto done; }
    for (size_t k = 1;; k <<= 1) {
#define KEY2(i) ((i) + k < n ? rk[(i) + k] + 1 : 0)
        size_t R = maxr + 2;
        memset(cnt, 0, R * sizeof *cnt);
        for (size_t i = 0; i < n; i++) cnt[KEY2(i)]++;
        for (size_t x = 1; x < R; x++) cnt[x] += cnt[x - 1];
        for (size_t i = n; i-- > 0;) tmp[--cnt[KEY2(i)]] = (uint32_t)i;
        memset(cnt, 0, (maxr + 1) * sizeof *cnt);
        for (size_t i = 0; i < n; i++) cnt[rk[i]]++;
        for (size_t x = 1; x <= maxr; x++) cnt[x] += cnt[x - 1];
        for (size_t q = n; q-- > 0;) { uint32_t i = tmp[q]; sa[--cnt[rk[i]]] = i; }
        tmp[sa[0]] = 0;
        for (size_t q = 1; q < n; q++) {
            uint32_t a = sa[q - 1], b = sa[q];
            tmp[b] = tmp[a] + !(rk[a] == rk[b] && KEY2(a) == KEY2(b));
        }
#undef KEY2
        uint32_t *s = rk; rk = tmp; tmp = s;
        maxr = rk[sa[n - 1]];
        if (maxr == n - 1) break;
    }
done:
    free(rk); free(tmp); free(cnt);
    return sa;
}

/* block RMQ over the LCP array (construction only) */
typedef struct {
    const uint32_t *a;
    size_t          n, nb;
    unsigned        J;
    uint32_t      **st;
} brmq_t;

static void brmq_init(brmq_t *q, const uint32_t *a, size_t n)
{
    q->a  = a; q->n = n; q->nb = (n + 63) / 64;
    q->J  = 1;
    while ((1ULL << q->J) <= q->nb) q->J++;
    q->st = malloc(q->J * sizeof *q->st);
    q->st[0] = malloc(q->nb * sizeof(uint32_t));
    for (size_t b = 0; b < q->nb; b++) {
        uint32_t m = UINT32_MAX;
        for (size_t i = b * 64; i < (b + 1) * 64 && i < n; i++) if (a[i] < m) m = a[i];
        q->st[0][b] = m;
    }
    for (unsigned j = 1; j < q->J; j++) {
        size_t w = (size_t)1 << j, cnt = q->nb >= w ? q->nb - w + 1 : 0;
        q->st[j] = malloc((cnt ? cnt : 1) * sizeof(uint32_t));
        for (size_t b = 0; b < cnt; b++) {
            uint32_t x = q->st[j - 1][b], y = q->st[j - 1][b + w / 2];
            q->st[j][b] = x < y ? x : y;
        }
    }
}

static void brmq_free(brmq_t *q)
{
    for (unsigned j = 0; j < q->J; j++) free(q->st[j]);
    free(q->st);
}

/* largest p <= x with a[p] < D (a[0] = 0 < D guarantees existence) */
static size_t left_less(const brmq_t *q, size_t x, uint32_t D)
{
    size_t bx = x / 64;
    for (size_t p = x;; p--) {
        if (q->a[p] < D) return p;
        if (p == bx * 64) break;
    }
    long cur = (long)bx - 1;
    for (unsigned j = q->J; j-- > 0;) {
        long w = 1L << j;
        if (cur + 1 >= w && q->st[j][cur + 1 - w] >= D) cur -= w;
    }
    size_t hi = (size_t)cur * 64 + 63;
    if (hi >= q->n) hi = q->n - 1;
    for (size_t p = hi;; p--) if (q->a[p] < D) return p;
}

/* smallest p >= s with a[p] < D, or n */
static size_t right_less(const brmq_t *q, size_t s, uint32_t D)
{
    if (s >= q->n) return q->n;
    size_t bs = s / 64, end = (bs + 1) * 64;
    for (size_t p = s; p < end && p < q->n; p++) if (q->a[p] < D) return p;
    size_t cur = bs + 1;
    for (unsigned j = q->J; j-- > 0;) {
        size_t w = (size_t)1 << j;
        if (cur + w <= q->nb && q->st[j][cur] >= D) cur += w;
    }
    if (cur >= q->nb) return q->n;
    for (size_t p = cur * 64; p < q->n; p++) if (q->a[p] < D) return p;
    return q->n;
}

typedef struct { uint32_t l, r, d, f; } snode_t;   /* f: 1 parens, 2 SDep */

static int cmp_open(const void *a, const void *b)
{
    const snode_t *x = a, *y = b;
    if (x->l != y->l) return x->l < y->l ? -1 : 1;
    if (x->r != y->r) return x->r > y->r ? -1 : 1;
    return 0;
}

static const snode_t *g_nodes;
static int cmp_close(const void *a, const void *b)
{
    const snode_t *x = &g_nodes[*(const uint32_t *)a], *y = &g_nodes[*(const uint32_t *)b];
    if (x->r != y->r) return x->r < y->r ? -1 : 1;
    if (x->l != y->l) return x->l > y->l ? -1 : 1;
    return 0;
}

static fcst_t *build(const uint8_t *text, size_t m, unsigned delta, unsigned delta_s,
                     unsigned sa_rate, int fast)
{
    if (m + 1 >= UINT32_MAX) return NULL;
    fcst_t *t = calloc(1, sizeof *t);
    size_t  n = m + 1;
    t->n       = n;
    t->h       = delta >= 2 ? delta / 2 : 1;
    t->hs      = delta_s >= 2 ? delta_s / 2 : 1;
    t->fast    = fast;
    t->sa_rate = sa_rate ? sa_rate : 32;

    int present[256] = { 0 };
    for (size_t i = 0; i < m; i++) {
        if (text[i] == 0) { free(t); return NULL; }
        present[text[i]] = 1;
    }
    unsigned k = 1;
    for (unsigned c = 1; c < 256; c++)
        if (present[c]) { t->map[c] = (uint8_t)k; t->sym[k] = (uint8_t)c; k++; }
    t->sigma = k;

    uint8_t *T = malloc(n);
    for (size_t i = 0; i < m; i++) T[i] = t->map[text[i]];
    T[m] = 0;

    uint32_t *SA = build_sa(T, n, t->sigma);

    /* C array and BWT */
    memset(t->C, 0, sizeof t->C);
    for (size_t i = 0; i < n; i++) t->C[T[i] + 1]++;
    for (unsigned c = 1; c <= t->sigma; c++) t->C[c] += t->C[c - 1];
    uint8_t *bw = malloc(n);
    for (size_t i = 0; i < n; i++) bw[i] = SA[i] ? T[SA[i] - 1] : T[n - 1];
    wm_build(&t->bwt, bw, n, t->sigma);
    free(bw);

    /* ISA, SA/ISA samples */
    uint32_t *ISA = malloc(n * sizeof *ISA);
    for (size_t i = 0; i < n; i++) ISA[SA[i]] = (uint32_t)i;
    bv_init(&t->sa_mark, n);
    size_t nmark = 0;
    for (size_t i = 0; i < n; i++)
        if (SA[i] % t->sa_rate == 0) { bv_set(&t->sa_mark, i); nmark++; }
    bv_build(&t->sa_mark);
    pa_init(&t->sa_samp, nmark, bits_needed(n / t->sa_rate));
    for (size_t i = 0, q = 0; i < n; i++)
        if (SA[i] % t->sa_rate == 0) pa_set(&t->sa_samp, q++, SA[i] / t->sa_rate);
    size_t nis = (n + t->sa_rate - 1) / t->sa_rate;
    pa_init(&t->isa_samp, nis, bits_needed(n));
    for (size_t q = 0; q < nis; q++) pa_set(&t->isa_samp, q, ISA[q * t->sa_rate]);

    /* LCP (Kasai) */
    uint32_t *LCP = malloc(n * sizeof *LCP);
    {
        size_t hh = 0;
        for (size_t i = 0; i < n; i++) {
            size_t r = ISA[i];
            if (r == 0) { LCP[0] = 0; hh = 0; continue; }
            size_t j = SA[r - 1];
            while (i + hh < n && j + hh < n && T[i + hh] == T[j + hh]) hh++;
            LCP[r] = (uint32_t)hh;
            if (hh) hh--;
        }
    }
    free(T);

    /* sampled nodes: SLink^h(v') for every internal v' with SDep = 0 mod h */
    brmq_t rq;
    brmq_init(&rq, LCP, n);
    size_t   cap = 1024, cnt = 0;
    snode_t *S = malloc(cap * sizeof *S);
    S[cnt++] = (snode_t){ 0, (uint32_t)(n - 1), 0, 3 };
    {
        uint32_t *stl = malloc((n + 1) * sizeof *stl), *stb = malloc((n + 1) * sizeof *stb);
        size_t    top = 0;
        stl[0] = 0; stb[0] = 0;
        for (size_t i = 1; i <= n; i++) {
            uint32_t cur = i < n ? LCP[i] : 0;
            uint32_t lb  = (uint32_t)(i - 1);
            while (cur < stl[top]) {
                uint32_t d = stl[top], l = stb[top], r = (uint32_t)(i - 1);
                top--;
                lb = l;
                /* SLink^hh(v') for hh = h (parentheses) and, in fast mode,
                   hh = hs (stored string depths) */
                for (int which = 0; which < (fast ? 2 : 1); which++) {
                    unsigned hh = which ? t->hs : t->h;
                    if (d % hh != 0 || d == hh) continue;
                    uint32_t D = d - hh;
                    size_t x = ISA[SA[l] + hh], y = ISA[SA[r] + hh];
                    size_t L = left_less(&rq, x, D);
                    size_t R = right_less(&rq, y + 1, D) - 1;
                    if (cnt == cap) { cap *= 2; S = realloc(S, cap * sizeof *S); }
                    S[cnt++] = (snode_t){ (uint32_t)L, (uint32_t)R, D,
                                          fast ? (which ? 2u : 1u) : 3u };
                }
            }
            if (cur > stl[top]) { ++top; stl[top] = cur; stb[top] = lb; }
        }
        free(stl); free(stb);
    }
    brmq_free(&rq);
    free(LCP); free(ISA); free(SA);

    qsort(S, cnt, sizeof *S, cmp_open);
    size_t u = 0;
    for (size_t i = 0; i < cnt; i++)
        if (u == 0 || S[i].l != S[u - 1].l || S[i].r != S[u - 1].r) S[u++] = S[i];
        else S[u - 1].f |= S[i].f;
    cnt = u;
    t->ns = cnt;

    if (!fast) {
        pa_init(&t->sd, cnt, bits_needed(n / t->h));
        for (size_t i = 0; i < cnt; i++) pa_set(&t->sd, i, S[i].d / t->h);
    } else {
        bv_init(&t->M, cnt);
        size_t bits = 0, k2 = 0;
        for (size_t i = 0; i < cnt; i++)
            if (S[i].f & 2) { bv_set(&t->M, i); bits += delta_len(S[i].d / t->hs + 1); k2++; }
        bv_build(&t->M);
        t->nsd    = k2;
        t->dsbits = bits;
        t->ds     = calloc(bits / 64 + 2, sizeof(uint64_t));
        pa_init(&t->dptr, k2 / DBLK + 1, bits_needed(bits));
        size_t pos = 0;
        k2 = 0;
        for (size_t i = 0; i < cnt; i++) {
            if (!(S[i].f & 2)) continue;
            if (k2 % DBLK == 0) pa_set(&t->dptr, k2 / DBLK, pos);
            put_delta(t->ds, &pos, S[i].d / t->hs + 1);
            k2++;
        }
    }

    /* interleave parentheses with leaves */
    uint32_t *cl = malloc(cnt * sizeof *cl);
    for (size_t i = 0; i < cnt; i++) cl[i] = (uint32_t)i;
    g_nodes = S;
    qsort(cl, cnt, sizeof *cl, cmp_close);
    t->np = 2 * cnt;
    bv_init(&t->B, n + t->np);
    bv_init(&t->P, t->np);
    {
        size_t bpos = 0, ppos = 0, oi = 0, ci = 0;
        for (size_t i = 0; i < n; i++) {
            while (oi < cnt && S[oi].l == i) { bv_set(&t->B, bpos++); bv_set(&t->P, ppos++); oi++; }
            bpos++;                                         /* leaf i */
            while (ci < cnt && S[cl[ci]].r == i) { bv_set(&t->B, bpos++); ppos++; ci++; }
        }
    }
    bv_build(&t->B);
    bv_build(&t->P);
    free(cl); free(S);

    /* block minima of the excess + segment tree */
    t->nb    = (t->np + 1 + 63) / 64;
    t->segsz = 1;
    while (t->segsz < t->nb) t->segsz <<= 1;
    t->seg = malloc(2 * t->segsz * sizeof *t->seg);
    for (size_t i = 0; i < 2 * t->segsz; i++) t->seg[i] = INT32_MAX;
    {
        long e = 0;
        for (size_t b = 0; b < t->nb; b++) {
            size_t k1 = b * 64 + 63;
            if (k1 > t->np) k1 = t->np;
            long mn = e;
            for (size_t k2 = b * 64; k2 < k1; k2++) {
                e += bv_get(&t->P, k2) ? 1 : -1;
                if (e < mn) mn = e;
            }
            if (k1 < t->np) e += bv_get(&t->P, k1) ? 1 : -1;   /* step into next block */
            t->seg[t->segsz + b] = (int32_t)mn;
        }
    }
    for (size_t i = t->segsz; i-- > 1;)
        t->seg[i] = t->seg[2 * i] < t->seg[2 * i + 1] ? t->seg[2 * i] : t->seg[2 * i + 1];
    return t;
}

fcst_t *fcst_build(const uint8_t *text, size_t m, unsigned delta, unsigned sa_rate)
{
    return build(text, m, delta, delta, sa_rate, 0);
}

fcst_t *fcst_build_fast(const uint8_t *text, size_t m, unsigned delta_par, unsigned delta_sdep,
                        unsigned sa_rate)
{
    return build(text, m, delta_par, delta_sdep, sa_rate, 1);
}

static size_t sdep_bytes(const fcst_t *t)
{
    if (!t->fast) return pa_bytes(&t->sd);
    return bv_bytes(&t->M) + (t->dsbits / 64 + 2) * 8 + pa_bytes(&t->dptr);
}

void fcst_free(fcst_t *t)
{
    if (!t) return;
    wm_free(&t->bwt);
    bv_free(&t->sa_mark); pa_free(&t->sa_samp); pa_free(&t->isa_samp);
    bv_free(&t->B); bv_free(&t->P); free(t->seg);
    if (t->fast) { bv_free(&t->M); free(t->ds); pa_free(&t->dptr); }
    else pa_free(&t->sd);
    free(t);
}

size_t fcst_bytes(const fcst_t *t)
{
    return sizeof *t + wm_bytes(&t->bwt) + bv_bytes(&t->sa_mark) + pa_bytes(&t->sa_samp) +
           pa_bytes(&t->isa_samp) + bv_bytes(&t->B) + bv_bytes(&t->P) +
           2 * t->segsz * sizeof(int32_t) + sdep_bytes(t);
}

void fcst_space(const fcst_t *t, FILE *o)
{
    double n = (double)(t->n);
#define ROW(name, bytes) fprintf(o, "  %-28s %12zu bytes  %7.3f bits/char\n", name, (size_t)(bytes), 8.0 * (bytes) / n)
    if (t->fast)
        fprintf(o, "FFCST: n = %zu, sigma = %u, delta_par = %u, delta_sdep = %u, sa_rate = %u, "
                   "sampled nodes = %zu, with SDep = %zu\n",
                t->n, t->sigma, 2 * t->h, 2 * t->hs, t->sa_rate, t->ns, t->nsd);
    else
        fprintf(o, "FCST: n = %zu, sigma = %u, delta = %u, sa_rate = %u, sampled nodes = %zu\n",
                t->n, t->sigma, 2 * t->h, t->sa_rate, t->ns);
    ROW("BWT wavelet matrix", wm_bytes(&t->bwt));
    ROW("SA samples (mark + values)", bv_bytes(&t->sa_mark) + pa_bytes(&t->sa_samp));
    ROW("ISA samples", pa_bytes(&t->isa_samp));
    ROW("sampled tree B + P", bv_bytes(&t->B) + bv_bytes(&t->P));
    ROW("excess min tree", 2 * t->segsz * sizeof(int32_t));
    if (t->fast) {
        ROW("SDep marks M", bv_bytes(&t->M));
        ROW("SDep Elias-delta + ptrs", (t->dsbits / 64 + 2) * 8 + pa_bytes(&t->dptr));
    } else
        ROW("sampled SDep / h", pa_bytes(&t->sd));
    ROW("total", fcst_bytes(t));
#undef ROW
}
