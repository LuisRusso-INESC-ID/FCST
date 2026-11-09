// FCST and Fast FCST on top of an SDSL compressed suffix array.
//
//   fcsts::fcst<t_csa> t;
//   t.build(file, delta_par, delta_sdep, fast, tmp_dir);   // original FCST: fast = false
//
// A node is its SA interval {l, r}. The interface follows SDSL's CSTs
// (root, parent, depth, sl, lca, select_child(v,1), sibling, child, lb, rb).
// The only CSA operations used are psi, lf (via backward_search),
// C / first-row symbol, csa[i] and csa.isa[i].
#pragma once

#include <sdsl/bit_vectors.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/construct_lcp.hpp>
#include <sdsl/suffix_arrays.hpp>

#include <algorithm>
#include <climits>
#include <string>
#include <vector>

namespace fcsts {

struct node {
    uint64_t l, r;
    bool operator==(const node &o) const { return l == o.l && r == o.r; }
    bool operator!=(const node &o) const { return !(*this == o); }
};

static const node NONE = { UINT64_MAX, UINT64_MAX };

template <class t_csa = sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<63>>, 32, 64>>
class fcst {
  public:
    typedef node     node_type;
    typedef uint64_t size_type;
    typedef typename t_csa::char_type char_type;

    t_csa csa;

  private:
    static const unsigned DBLK = 16;     // Elias-delta pointer sampling
    static const unsigned DTHR = 64;     // child(): backward mapping if SDep <= DTHR

    uint64_t n = 0;
    bool     fast = false;
    unsigned h = 1, hs = 1;              // delta_par / 2, delta_sdep / 2
    uint64_t ns = 0, nsd = 0, np = 0, nb = 0, segsz = 1;

    sdsl::bit_vector                 B, P, M;
    sdsl::rank_support_v5<1>         P_rank, M_rank;
    sdsl::select_support_mcl<0>      B_sel0;
    sdsl::select_support_mcl<1>      B_sel1;
    sdsl::int_vector<>               seg;      // min excess per block of 64 (segment tree)
    uint64_t                         seg_inf = 0;
    sdsl::int_vector<>               sd;       // original: SDep / h per sampled node
    std::vector<uint64_t>            ds;       // fast: Elias-delta codes of SDep / hs + 1
    uint64_t                         dsbits = 0;
    sdsl::int_vector<>               dptr;

    fcst(const fcst &) = delete;
    fcst &operator=(const fcst &) = delete;

  public:
    fcst() {}

    // ------------------------------------------------------------ queries
    size_type size() const { return n; }
    node root() const { return node{ 0, n - 1 }; }
    node select_leaf(size_type i) const { return node{ i - 1, i - 1 }; }   // 1-based, as in SDSL
    size_type lb(node v) const { return v.l; }
    size_type rb(node v) const { return v.r; }
    static bool is_leaf(node v) { return v.l == v.r; }
    static bool ancestor(node v, node w) { return v.l <= w.l && w.r <= v.r; }
    uint64_t sampled() const { return ns; }
    uint64_t sampled_sdep() const { return fast ? nsd : ns; }

    size_type depth(node v) const
    {
        if (is_leaf(v)) return n - csa[v.l];
        return fast ? sdep_fast(v.l, v.r) : lca_orig(v.l, v.r, nullptr);
    }

    node lca(node v, node w) const
    {
        if (ancestor(v, w)) return v;
        if (ancestor(w, v)) return w;
        return lca_node(std::min(v.l, w.l), std::max(v.r, w.r));
    }

    node parent(node v) const
    {
        if (v == root()) return root();          // SDSL convention
        node best = NONE;
        if (v.l > 0) best = lca_node(v.l - 1, v.l);
        if (v.r < n - 1) {
            node y = lca_node(v.r, v.r + 1);
            if (best == NONE || y.r - y.l < best.r - best.l) best = y;
        }
        return best;
    }

    node sl(node v) const
    {
        if (v == root()) return root();
        if (is_leaf(v)) {
            if (v.l == 0) return root();
            uint64_t x = csa.psi[v.l];
            return node{ x, x };
        }
        return lca_node(csa.psi[v.l], csa.psi[v.r]);
    }

    node select_child(node v, size_type) const   // first child
    {
        if (is_leaf(v)) return root();
        return child_range(v.l, v.r, depth(v), -1);
    }

    node sibling(node v) const
    {
        if (v == root()) return root();
        node p = parent(v);
        if (v.r == p.r) return root();
        return child_range(v.r + 1, p.r, depth(p), -1);
    }

    node child(node v, char_type c) const
    {
        if (is_leaf(v)) return root();
        uint64_t cc = csa.char2comp[c];
        if (cc == 0 && c != 0) return root();
        node w = child_range(v.l, v.r, depth(v), (int64_t)cc);
        return w == NONE ? root() : w;
    }

    // ------------------------------------------------------------ space
    uint64_t bytes() const
    {
        return sdsl::size_in_bytes(csa) + tree_bytes();
    }
    uint64_t tree_bytes() const
    {
        uint64_t t = sdsl::size_in_bytes(B) + sdsl::size_in_bytes(P) + sdsl::size_in_bytes(P_rank) +
                     sdsl::size_in_bytes(B_sel0) + sdsl::size_in_bytes(B_sel1) + sdsl::size_in_bytes(seg);
        if (fast)
            t += sdsl::size_in_bytes(M) + sdsl::size_in_bytes(M_rank) + ds.size() * 8 + sdsl::size_in_bytes(dptr);
        else
            t += sdsl::size_in_bytes(sd);
        return t;
    }

  private:
    // ------------------------------------------------------------ CSA helpers
    uint64_t fcode(uint64_t i) const     // comp. symbol of F[i]
    {
        uint64_t lo = 0, hi = csa.sigma - 1;
        while (lo < hi) {
            uint64_t mid = (lo + hi + 1) / 2;
            if (csa.C[mid] <= i) lo = mid; else hi = mid - 1;
        }
        return lo;
    }

    bool lf_interval(uint64_t comp, uint64_t &l, uint64_t &r) const
    {
        uint64_t nl, nr;
        if (!sdsl::backward_search(csa, l, r, csa.comp2char[comp], nl, nr)) return false;
        l = nl; r = nr;
        return true;
    }

    uint64_t psi_pow(uint64_t x, uint64_t d) const   // Psi^d(x)
    {
        if (d <= 32) { while (d--) x = csa.psi[x]; return x; }
        return csa.isa[csa[x] + d];
    }

    // ------------------------------------------------------------ parentheses
    long excess(uint64_t k) const { return 2 * (long)P_rank(k) - (long)k; }

    long scan_min(uint64_t k0, uint64_t k1) const
    {
        long e = excess(k0), m = e;
        for (uint64_t k = k0; k < k1; k++) { e += P[k] ? 1 : -1; if (e < m) m = e; }
        return m;
    }

    long seg_min(uint64_t L, uint64_t R) const
    {
        uint64_t m = seg_inf;
        L += segsz; R += segsz + 1;
        while (L < R) {
            if (L & 1) m = std::min<uint64_t>(m, seg[L++]);
            if (R & 1) m = std::min<uint64_t>(m, seg[--R]);
            L >>= 1; R >>= 1;
        }
        return (long)m;
    }

    long seg_right(uint64_t nd, uint64_t nl, uint64_t nr, uint64_t R, long th) const
    {
        if (nl >= R || (long)seg[nd] > th) return -1;
        if (nr - nl == 1) return (long)nl;
        uint64_t mid = (nl + nr) / 2;
        long     x   = seg_right(2 * nd + 1, mid, nr, R, th);
        return x >= 0 ? x : seg_right(2 * nd, nl, mid, R, th);
    }

    long seg_left(uint64_t nd, uint64_t nl, uint64_t nr, uint64_t L, long th) const
    {
        if (nr <= L || (long)seg[nd] > th) return -1;
        if (nr - nl == 1) return (long)nl;
        uint64_t mid = (nl + nr) / 2;
        long     x   = seg_left(2 * nd, nl, mid, L, th);
        return x >= 0 ? x : seg_left(2 * nd + 1, mid, nr, L, th);
    }

    long rmq_min(uint64_t a, uint64_t b) const
    {
        uint64_t ba = a / 64, bb = b / 64;
        if (ba == bb) return scan_min(a, b);
        long m = scan_min(a, ba * 64 + 63);
        if (ba + 1 <= bb - 1) m = std::min(m, seg_min(ba + 1, bb - 1));
        return std::min(m, scan_min(bb * 64, b));
    }

    uint64_t bwd_search(uint64_t a, long th) const      // largest p < a, E(p) <= th
    {
        uint64_t p = a, start = (a / 64) * 64;
        long     e = excess(a);
        while (p > start) { p--; e -= P[p] ? 1 : -1; if (e <= th) return p; }
        long     b  = seg_right(1, 0, segsz, a / 64, th);
        uint64_t k0 = (uint64_t)b * 64, k1 = std::min<uint64_t>(k0 + 63, np), best = k0;
        e = excess(k0);
        for (uint64_t k = k0; k < k1; k++) { e += P[k] ? 1 : -1; if (e <= th) best = k + 1; }
        return best;
    }

    uint64_t fwd_search(uint64_t p, long th) const      // smallest k > p, E(k) <= th
    {
        uint64_t k = p, end = std::min<uint64_t>((p / 64) * 64 + 63, np);
        long     e = excess(p);
        while (k < end) { e += P[k] ? 1 : -1; k++; if (e <= th) return k; }
        long b = seg_left(1, 0, segsz, p / 64 + 1, th);
        k = (uint64_t)b * 64;
        e = excess(k);
        while (e > th) { e += P[k] ? 1 : -1; k++; }
        return k;
    }

    void sampled_interval(uint64_t p, uint64_t &l, uint64_t &r) const
    {
        uint64_t q = fwd_search(p, excess(p)) - 1;
        l = B_sel1(p + 1) - p;
        r = B_sel1(q + 1) - q - 1;
    }

    uint64_t lcsa(uint64_t i, uint64_t j) const         // '(' position in P
    {
        uint64_t a = B_sel0(i + 1) - i, b = B_sel0(j + 1) - j;
        return bwd_search(a, rmq_min(a, b) - 1);
    }

    // ------------------------------------------------------------ SDep storage
    static int  gbit(const std::vector<uint64_t> &w, uint64_t p) { return (w[p >> 6] >> (p & 63)) & 1; }
    static void pbits(std::vector<uint64_t> &w, uint64_t &p, uint64_t v, unsigned nbits)
    {
        for (unsigned q = nbits; q-- > 0; p++) if ((v >> q) & 1) w[p >> 6] |= 1ULL << (p & 63);
    }
    static unsigned dlen(uint64_t x)
    {
        unsigned N = 63 - __builtin_clzll(x), L = N + 1, LL = 31 - __builtin_clz(L);
        return 2 * LL + 1 + N;
    }
    static void pdelta(std::vector<uint64_t> &w, uint64_t &p, uint64_t x)
    {
        unsigned N = 63 - __builtin_clzll(x), L = N + 1, LL = 31 - __builtin_clz(L);
        p += LL;
        pbits(w, p, L, LL + 1);
        pbits(w, p, x, N);
    }
    uint64_t gdelta(uint64_t &p) const
    {
        unsigned LL = 0;
        while (!gbit(ds, p)) { LL++; p++; }
        uint64_t L = 0, x = 1;
        for (unsigned q = 0; q <= LL; q++) L = (L << 1) | gbit(ds, p++);
        for (uint64_t q = 1; q < L; q++) x = (x << 1) | gbit(ds, p++);
        return x;
    }
    long stored_sdep(uint64_t id) const
    {
        if (!fast) return (long)sd[id] * h;
        if (!M[id]) return -1;
        uint64_t rk = M_rank(id), pos = dptr[rk / DBLK];
        for (uint64_t q = rk % DBLK; q > 0; q--) gdelta(pos);
        return (long)(gdelta(pos) - 1) * hs;
    }

    // ------------------------------------------------------------ LCA / SDep
    uint64_t lca_orig(uint64_t i, uint64_t j, node *out) const   // original FCST
    {
        unsigned lim = 2 * h;
        uint64_t let[512];
        std::vector<uint64_t> big;
        uint64_t *L = lim <= 512 ? let : (big.resize(lim), big.data());
        uint64_t  best = 0, bestk = 0, bestp = 0;
        for (unsigned k = 0; k < lim; k++) {
            uint64_t ci = fcode(i), cj = fcode(j);
            if (ci != cj) { if (k >= best) { best = k; bestk = k; bestp = 0; } break; }
            uint64_t p = lcsa(i, j), cand = k + (uint64_t)stored_sdep(P_rank(p));
            if (cand > best || k == 0) { best = cand; bestk = k; bestp = p; }
            L[k] = ci;
            i = csa.psi[i]; j = csa.psi[j];
        }
        if (out) {
            uint64_t l, r;
            sampled_interval(bestp, l, r);
            for (uint64_t q = bestk; q-- > 0;) lf_interval(L[q], l, r);
            *out = node{ l, r };
        }
        return best;
    }

    node lca_fast(uint64_t i, uint64_t j) const   // smallest LF-mapped LCSA
    {
        unsigned lim = 2 * h;
        uint64_t lb_[256], pb_[256];
        std::vector<uint64_t> bl, bp;
        uint64_t *let = lim <= 256 ? lb_ : (bl.resize(lim), bl.data());
        uint64_t *par = lim <= 256 ? pb_ : (bp.resize(lim), bp.data());
        unsigned kk = 0;
        for (unsigned k = 0; k < lim; k++) {
            uint64_t ci = fcode(i), cj = fcode(j);
            kk = k + 1;
            if (ci != cj) { par[k] = 0; break; }
            par[k] = lcsa(i, j);
            let[k] = ci;
            i = csa.psi[i]; j = csa.psi[j];
        }
        uint64_t l, r, l2, r2;
        sampled_interval(par[kk - 1], l, r);
        for (unsigned k = kk - 1; k-- > 0;) {
            lf_interval(let[k], l, r);
            if (par[k] == 0) continue;
            sampled_interval(par[k], l2, r2);
            if (r2 - l2 < r - l) { l = l2; r = r2; }
        }
        return node{ l, r };
    }

    uint64_t sdep_fast(uint64_t i, uint64_t j) const
    {
        unsigned lim = 2 * hs;
        uint64_t best = 0;
        for (unsigned k = 0; k < lim; k++) {
            if (fcode(i) != fcode(j)) { best = std::max<uint64_t>(best, k); break; }
            long s = stored_sdep(P_rank(lcsa(i, j)));
            if (s >= 0) best = std::max<uint64_t>(best, k + (uint64_t)s);
            i = csa.psi[i]; j = csa.psi[j];
        }
        return best;
    }

    node lca_node(uint64_t i, uint64_t j) const
    {
        if (fast) return lca_fast(i, j);
        node v;
        lca_orig(i, j, &v);
        return v;
    }

    // Children of the node whose leaves lo..hi share a prefix of length d.
    // c >= 0: child by comp. symbol c; c < 0: child starting at leaf lo.
    node child_range(uint64_t lo, uint64_t hi, uint64_t d, int64_t c) const
    {
        if (d <= DTHR) {   // Psi^d keeps the order; map the c-range back with LF
            uint64_t let[DTHR];
            uint64_t a = lo, b = hi;
            for (uint64_t k = 0; k < d; k++) { let[k] = fcode(a); a = csa.psi[a]; b = csa.psi[b]; }
            uint64_t cc = c >= 0 ? (uint64_t)c : fcode(a);
            uint64_t l = std::max<uint64_t>(a, csa.C[cc]), r = std::min<uint64_t>(b, csa.C[cc + 1] - 1);
            if (l > r) return NONE;
            for (uint64_t k = d; k-- > 0;) if (!lf_interval(let[k], l, r)) return NONE;
            return node{ l, r };
        }
        uint64_t cc = c >= 0 ? (uint64_t)c : fcode(psi_pow(lo, d));
        uint64_t x = lo, y = hi + 1;
        while (x < y) { uint64_t m = x + (y - x) / 2; if (fcode(psi_pow(m, d)) < cc) x = m + 1; else y = m; }
        uint64_t s = x;
        y = hi + 1;
        while (x < y) { uint64_t m = x + (y - x) / 2; if (fcode(psi_pow(m, d)) <= cc) x = m + 1; else y = m; }
        if (s >= x) return NONE;
        return node{ s, x - 1 };
    }

    // ------------------------------------------------------------ construction
    struct brmq {
        const sdsl::int_vector<> *a;
        uint64_t n, nb;
        unsigned J;
        std::vector<std::vector<uint32_t>> st;
        brmq(const sdsl::int_vector<> &v) : a(&v), n(v.size()), nb((v.size() + 63) / 64)
        {
            J = 1;
            while ((1ULL << J) <= nb) J++;
            st.resize(J);
            st[0].resize(nb);
            for (uint64_t b = 0; b < nb; b++) {
                uint64_t m = UINT64_MAX;
                for (uint64_t i = b * 64; i < (b + 1) * 64 && i < n; i++) m = std::min<uint64_t>(m, v[i]);
                st[0][b] = (uint32_t)std::min<uint64_t>(m, UINT32_MAX);
            }
            for (unsigned j = 1; j < J; j++) {
                uint64_t w = 1ULL << j;
                st[j].resize(nb >= w ? nb - w + 1 : 1);
                for (uint64_t b = 0; b + w <= nb; b++) st[j][b] = std::min(st[j - 1][b], st[j - 1][b + w / 2]);
            }
        }
        uint64_t left_less(uint64_t x, uint64_t D) const
        {
            uint64_t bx = x / 64;
            for (uint64_t p = x;; p--) { if ((*a)[p] < D) return p; if (p == bx * 64) break; }
            long cur = (long)bx - 1;
            for (unsigned j = J; j-- > 0;) {
                long w = 1L << j;
                if (cur + 1 >= w && st[j][cur + 1 - w] >= D) cur -= w;
            }
            for (uint64_t p = std::min<uint64_t>((uint64_t)cur * 64 + 63, n - 1);; p--) if ((*a)[p] < D) return p;
        }
        uint64_t right_less(uint64_t s, uint64_t D) const
        {
            if (s >= n) return n;
            uint64_t bs = s / 64;
            for (uint64_t p = s; p < (bs + 1) * 64 && p < n; p++) if ((*a)[p] < D) return p;
            uint64_t cur = bs + 1;
            for (unsigned j = J; j-- > 0;) {
                uint64_t w = 1ULL << j;
                if (cur + w <= nb && st[j][cur] >= D) cur += w;
            }
            if (cur >= nb) return n;
            for (uint64_t p = cur * 64; p < n; p++) if ((*a)[p] < D) return p;
            return n;
        }
    };

    struct snode { uint64_t l, r, d; unsigned f; };

  public:
    void build(const std::string &file, unsigned delta_par, unsigned delta_sdep, bool fast_mode,
               const std::string &tmp_dir = "/tmp")
    {
        fast = fast_mode;
        h    = delta_par >= 2 ? delta_par / 2 : 1;
        hs   = delta_sdep >= 2 ? delta_sdep / 2 : 1;
        if (!fast) hs = h;

        sdsl::cache_config cfg(false, tmp_dir, "fcst_" + std::to_string(getpid()));
        sdsl::construct(csa, file, cfg, 1);
        sdsl::construct_lcp_kasai<8>(cfg);
        n = csa.size();

        std::vector<snode> S;
        {
            sdsl::int_vector<> SA, LCP;
            sdsl::load_from_cache(SA, sdsl::conf::KEY_SA, cfg);
            sdsl::load_from_cache(LCP, sdsl::conf::KEY_LCP, cfg);
            sdsl::int_vector<> ISA(n, 0, sdsl::bits::hi(n) + 1);
            for (uint64_t i = 0; i < n; i++) ISA[SA[i]] = i;
            brmq rq(LCP);
            S.push_back(snode{ 0, n - 1, 0, 3 });
            std::vector<uint64_t> stl{ 0 }, stb{ 0 };
            for (uint64_t i = 1; i <= n; i++) {
                uint64_t cur = i < n ? (uint64_t)LCP[i] : 0, lb = i - 1;
                while (cur < stl.back()) {
                    uint64_t d = stl.back(), l = stb.back(), r = i - 1;
                    stl.pop_back(); stb.pop_back();
                    lb = l;
                    for (int which = 0; which < (fast ? 2 : 1); which++) {
                        unsigned hh = which ? hs : h;
                        if (d % hh != 0 || d == hh) continue;
                        uint64_t D = d - hh;
                        uint64_t x = ISA[SA[l] + hh], y = ISA[SA[r] + hh];
                        uint64_t L = rq.left_less(x, D), R = rq.right_less(y + 1, D) - 1;
                        S.push_back(snode{ L, R, D, fast ? (which ? 2u : 1u) : 3u });
                    }
                }
                if (cur > stl.back()) { stl.push_back(cur); stb.push_back(lb); }
            }
        }
        sdsl::util::delete_all_files(cfg.file_map);

        std::sort(S.begin(), S.end(), [](const snode &a, const snode &b) {
            return a.l != b.l ? a.l < b.l : a.r > b.r;
        });
        uint64_t u = 0;
        for (uint64_t i = 0; i < S.size(); i++) {
            if (u == 0 || S[i].l != S[u - 1].l || S[i].r != S[u - 1].r) S[u++] = S[i];
            else S[u - 1].f |= S[i].f;
        }
        S.resize(u);
        ns = u;

        if (!fast) {
            sd = sdsl::int_vector<>(ns, 0, sdsl::bits::hi(n / h) + 1);
            for (uint64_t i = 0; i < ns; i++) sd[i] = S[i].d / h;
        } else {
            M = sdsl::bit_vector(ns, 0);
            uint64_t bits = 0;
            nsd = 0;
            for (uint64_t i = 0; i < ns; i++)
                if (S[i].f & 2) { M[i] = 1; bits += dlen(S[i].d / hs + 1); nsd++; }
            sdsl::util::init_support(M_rank, &M);
            dsbits = bits;
            ds.assign(bits / 64 + 2, 0);
            dptr = sdsl::int_vector<>(nsd / DBLK + 1, 0, sdsl::bits::hi(bits + 1) + 1);
            uint64_t pos = 0, k = 0;
            for (uint64_t i = 0; i < ns; i++) {
                if (!(S[i].f & 2)) continue;
                if (k % DBLK == 0) dptr[k / DBLK] = pos;
                pdelta(ds, pos, S[i].d / hs + 1);
                k++;
            }
        }

        std::vector<uint32_t> cl(ns);
        for (uint64_t i = 0; i < ns; i++) cl[i] = (uint32_t)i;
        std::sort(cl.begin(), cl.end(), [&](uint32_t a, uint32_t b) {
            return S[a].r != S[b].r ? S[a].r < S[b].r : S[a].l > S[b].l;
        });
        np = 2 * ns;
        B  = sdsl::bit_vector(n + np, 0);
        P  = sdsl::bit_vector(np, 0);
        {
            uint64_t bpos = 0, ppos = 0, oi = 0, ci = 0;
            for (uint64_t i = 0; i < n; i++) {
                while (oi < ns && S[oi].l == i) { B[bpos++] = 1; P[ppos++] = 1; oi++; }
                bpos++;
                while (ci < ns && S[cl[ci]].r == i) { B[bpos++] = 1; ppos++; ci++; }
            }
        }
        sdsl::util::init_support(P_rank, &P);
        sdsl::util::init_support(B_sel0, &B);
        sdsl::util::init_support(B_sel1, &B);

        nb    = (np + 1 + 63) / 64;
        segsz = 1;
        while (segsz < nb) segsz <<= 1;
        unsigned w = sdsl::bits::hi(ns + 2) + 2;
        seg_inf    = (1ULL << w) - 1;
        seg        = sdsl::int_vector<>(2 * segsz, seg_inf, w);
        long e = 0;
        for (uint64_t b = 0; b < nb; b++) {
            uint64_t k1 = std::min<uint64_t>(b * 64 + 63, np);
            long     mn = e;
            for (uint64_t k = b * 64; k < k1; k++) { e += P[k] ? 1 : -1; mn = std::min(mn, e); }
            if (k1 < np) e += P[k1] ? 1 : -1;
            seg[segsz + b] = (uint64_t)mn;
        }
        for (uint64_t i = segsz; i-- > 1;) seg[i] = std::min<uint64_t>(seg[2 * i], seg[2 * i + 1]);
    }
};

}   // namespace fcsts
