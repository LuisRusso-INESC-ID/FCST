// Unified benchmark for compressed suffix trees.
// usage: bench_cst FILE STRUCT [verify]
//   STRUCT: sada | sct3 | sct3-t2 | cn | cnr | fcst:D | ffcst:DP:DS
// Methodology (Gog's thesis / Canovas-Navarro):
//   V1: ancestors of the parents of random leaves  -> parent, depth, fchild, sibling, child
//   V2: suffix-link chains from parents of random leaves -> sl
//   V3: random leaf pairs -> lca
#include <sys/resource.h>

#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sdsl/suffix_trees.hpp>
#include <string>
#include <vector>

#include "cst_cn.h"
#include "fcst_sdsl.hpp"
extern "C" {
#include "fcst.h"
}

using namespace std;
using timer = chrono::high_resolution_clock;

static const size_t V1_LEAVES = 10000, V2_LEAVES = 1000, V3_PAIRS = 100000, CAP = 100000;

struct rng_t {
    uint64_t s = 0x2545F4914F6CDD1DULL;
    uint64_t operator()() { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return s; }
};

static string text;

// ---- FCST adapter with the sdsl interface subset used here ------------------
struct fcst_cst {
    typedef fcst_node node_type;
    fcst_t *t = nullptr;
    fcst_node root() const { return fcst_root(t); }
    size_t size() const { return fcst_n(t); }
    fcst_node select_leaf(size_t i) const { fcst_node v = { i - 1, i - 1 }; return v; }
    fcst_node parent(fcst_node v) const { fcst_node p = fcst_parent(t, v); return fcst_is_none(p) ? root() : p; }
    size_t depth(fcst_node v) const { return fcst_sdep(t, v); }
    fcst_node sl(fcst_node v) const { return fcst_slink(t, v); }
    fcst_node lca(fcst_node a, fcst_node b) const { return fcst_lca(t, a, b); }
    fcst_node nz(fcst_node v) const { return fcst_is_none(v) ? root() : v; }
    fcst_node select_child(fcst_node v, size_t) const { return nz(fcst_fchild(t, v)); }
    fcst_node sibling(fcst_node v) const { return nz(fcst_nsibling(t, v)); }
    fcst_node child(fcst_node v, unsigned char c) const { return nz(fcst_child(t, v, c)); }
    size_t lb(fcst_node v) const { return v.l; }
    size_t rb(fcst_node v) const { return v.r; }
};
static bool operator==(fcst_node a, fcst_node b) { return a.l == b.l && a.r == b.r; }
static bool operator!=(fcst_node a, fcst_node b) { return !(a == b); }
static size_t bytes_of(const fcst_cst &c) { return fcst_bytes(c.t); }
template <class C> static size_t bytes_of(const C &c) { return sdsl::size_in_bytes(c); }
template <class A> static size_t bytes_of(const fcsts::fcst<A> &c) { return c.bytes(); }
typedef sdsl::csa_wt<sdsl::wt_huff<sdsl::rrr_vector<63>>, 32, 64> csa_rrr;

static double peak_mb()
{
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
    return ru.ru_maxrss / 1024.0;
}

template <class F> static double ns_per(size_t cnt, F f)
{
    auto a = timer::now();
    f();
    auto b = timer::now();
    return cnt ? chrono::duration_cast<chrono::nanoseconds>(b - a).count() / (double)cnt : 0;
}

template <class C> static void run(const string &name, C &cst, double build_s, double build_mb, bool verify)
{
    typedef typename C::node_type node_t;
    size_t n = cst.size();
    auto   root = cst.root();
    volatile size_t sink = 0;

    // V1
    rng_t         r1;
    vector<node_t> V1;
    vector<unsigned char> sym;
    for (size_t q = 0; q < V1_LEAVES && V1.size() < CAP; q++) {
        auto v = cst.parent(cst.select_leaf(1 + r1() % (n - 1)));
        while (v != root && V1.size() < CAP) {
            V1.push_back(v);
            sym.push_back((unsigned char)text[r1() % text.size()]);
            v = cst.parent(v);
        }
    }
    double t_par = ns_per(V1.size(), [&] { for (auto &v : V1) sink += cst.lb(cst.parent(v)); });
    double t_dep = ns_per(V1.size(), [&] { for (auto &v : V1) sink += cst.depth(v); });
    double t_fch = ns_per(V1.size(), [&] { for (auto &v : V1) sink += cst.lb(cst.select_child(v, 1)); });
    double t_sib = ns_per(V1.size(), [&] { for (auto &v : V1) sink += cst.lb(cst.sibling(v)); });
    double t_chi = ns_per(V1.size(), [&] { for (size_t i = 0; i < V1.size(); i++) sink += cst.lb(cst.child(V1[i], sym[i])); });

    // V2
    rng_t         r2; r2.s ^= 0x1234;
    vector<node_t> V2;
    for (size_t q = 0; q < V2_LEAVES && V2.size() < CAP; q++) {
        auto v = cst.parent(cst.select_leaf(1 + r2() % (n - 1)));
        while (v != root && V2.size() < CAP) { V2.push_back(v); v = cst.sl(v); }
    }
    double t_sl = ns_per(V2.size(), [&] { for (auto &v : V2) sink += cst.lb(cst.sl(v)); });

    // V3
    rng_t         r3; r3.s ^= 0x9876;
    vector<node_t> V3;
    for (size_t q = 0; q < V3_PAIRS; q++) {
        size_t a = 1 + r3() % (n - 1), b = 1 + r3() % (n - 1);
        if (a > b) swap(a, b);
        V3.push_back(cst.select_leaf(a));
        V3.push_back(cst.select_leaf(b));
    }
    double t_lca = ns_per(V3_PAIRS, [&] { for (size_t i = 0; i < V3.size(); i += 2) sink += cst.lb(cst.lca(V3[i], V3[i + 1])); });

    size_t bytes = bytes_of(cst);
    printf("RESULT\t%s\t%zu\t%.3f\t%.2f\t%.1f\t%zu\t%zu\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n", name.c_str(), n,
           bytes * 8.0 / n, build_s, build_mb, V1.size(), V2.size(), t_par, t_dep, t_fch, t_sib, t_chi, t_sl, t_lca);

    if (verify) {   // dump (lb, rb, depth) of the first V1 nodes for cross-checking
        FILE *f = fopen((name + ".verify").c_str(), "w");
        for (size_t i = 0; i < V1.size() && i < 2000; i++)
            fprintf(f, "%zu %zu %zu\n", (size_t)cst.lb(V1[i]), (size_t)cst.rb(V1[i]), (size_t)cst.depth(V1[i]));
        for (size_t i = 0; i < V1.size() && i < 2000; i++) {
            auto a = cst.select_child(V1[i], 1), b = cst.sibling(V1[i]), c = cst.child(V1[i], sym[i]);
            fprintf(f, "nav %zu %zu %zu %zu %zu %zu\n", (size_t)cst.lb(a), (size_t)cst.rb(a), (size_t)cst.lb(b),
                    (size_t)cst.rb(b), (size_t)cst.lb(c), (size_t)cst.rb(c));
        }
        for (size_t i = 0; i < V2.size() && i < 2000; i++)
            fprintf(f, "sl %zu %zu\n", (size_t)cst.lb(V2[i]), (size_t)cst.rb(V2[i]));
        for (size_t i = 0; i < 4000; i += 2) {
            auto w = cst.lca(V3[i], V3[i + 1]);
            fprintf(f, "lca %zu %zu\n", (size_t)cst.lb(w), (size_t)cst.rb(w));
        }
        fclose(f);
    }
    (void)sink;
}

template <class C> static void run_sdsl(const string &name, const string &file, bool verify)
{
    C                  cst;
    sdsl::cache_config cfg(true, "/tmp/sdsltmp", name + "_" + to_string(getpid()));
    auto a = timer::now();
    sdsl::construct(cst, file, cfg, 1);
    double bs = chrono::duration<double>(timer::now() - a).count();
    run(name, cst, bs, peak_mb(), verify);
}

int main(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: %s FILE STRUCT [verify]\n", argv[0]); return 1; }
    string file = argv[1], s = argv[2];
    bool   verify = argc > 3;
    {
        ifstream in(file, ios::binary);
        text.assign(istreambuf_iterator<char>(in), istreambuf_iterator<char>());
    }
    using namespace sdsl;
    if (s == "sada")         run_sdsl<cst_sada<>>(s, file, verify);
    else if (s == "sct3")    run_sdsl<cst_sct3<>>(s, file, verify);
    else if (s == "sct3-t2") run_sdsl<cst_sct3<csa_wt<>, lcp_support_tree2<>>>(s, file, verify);
    else if (s == "cn")      run_sdsl<cstds::cst_cn<csa_wt<>, lcp_dac<>, cstds::npr_support_cn, 32>>(s, file, verify);
    else if (s == "cnr")     run_sdsl<cstds::cst_cn<csa_wt<>, lcp_dac<>, cstds::npr_support_cnr, 32, 8>>(s, file, verify);
    else if (s == "sada-rrr")    run_sdsl<cst_sada<csa_rrr>>(s, file, verify);
    else if (s == "sct3-rrr")    run_sdsl<cst_sct3<csa_rrr>>(s, file, verify);
    else if (s == "sct3-t2-rrr") run_sdsl<cst_sct3<csa_rrr, lcp_support_tree2<>>>(s, file, verify);
    else if (s == "cn-rrr")      run_sdsl<cstds::cst_cn<csa_rrr, lcp_dac<>, cstds::npr_support_cn, 32>>(s, file, verify);
    else if (s == "cnr-rrr")     run_sdsl<cstds::cst_cn<csa_rrr, lcp_dac<>, cstds::npr_support_cnr, 32, 8>>(s, file, verify);
    else if (s.rfind("fcst-rrr:", 0) == 0 || s.rfind("ffcst-rrr:", 0) == 0) {
        unsigned dp = 0, ds = 0;
        bool     fast = s[1] == 'f';
        if (fast) sscanf(s.c_str(), "ffcst-rrr:%u:%u", &dp, &ds);
        else      sscanf(s.c_str(), "fcst-rrr:%u", &dp);
        fcsts::fcst<csa_rrr> c;
        auto a = timer::now();
        c.build(file, dp, fast ? ds : dp, fast, "/tmp/sdsltmp");
        double bs = chrono::duration<double>(timer::now() - a).count();
        run(s, c, bs, peak_mb(), verify);
    }
    else if (s.rfind("fcst:", 0) == 0 || s.rfind("ffcst:", 0) == 0) {
        fcst_cst c;
        unsigned dp = 0, ds = 0;
        bool     fast = s[0] == 'f' && s[1] == 'f';
        if (fast) sscanf(s.c_str(), "ffcst:%u:%u", &dp, &ds);
        else      sscanf(s.c_str(), "fcst:%u", &dp);
        auto a = timer::now();
        const uint8_t *p = (const uint8_t *)text.data();
        c.t = fast ? fcst_build_fast(p, text.size(), dp, ds, 32) : fcst_build(p, text.size(), dp, 32);
        double bs = chrono::duration<double>(timer::now() - a).count();
        run(s, c, bs, peak_mb(), verify);
        fcst_free(c.t);
    } else { fprintf(stderr, "unknown structure %s\n", s.c_str()); return 1; }
    return 0;
}
