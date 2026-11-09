#!/usr/bin/env python3
"""Regenerate tab_same.tex, tab_pizza.tex, tab_default.tex, tab_build.tex and p_*.dat from bench_cst TSV output.
Texts in PIZZA go to tab_pizza.tex; all others (the 8 MB texts) to tab_same.tex. Tables hold at most
PER_FLOAT texts per float; further floats get the labels <label>-2, <label>-3, ...
usage: gen_tables.py [results_same_csa.tsv results_default.tsv [OUTDIR]]
Defaults: results/results_same_csa.tsv results/results_default.tsv ../paper/generated
(paths relative to the repository root)."""
import os, sys, collections
cols = ["text","s","n","bpc","bs","mb","v1","v2","par","dep","fch","sib","chi","sl","lca"]
def load(fn):
    return [dict(zip(cols, l.rstrip("\n").split("\t"))) for l in open(fn) if not l.startswith("text")]
here = os.path.dirname(os.path.abspath(__file__))
code = os.path.dirname(here)
a_same = sys.argv[1] if len(sys.argv) > 1 else os.path.join(code, "results", "results_same_csa.tsv")
a_dflt = sys.argv[2] if len(sys.argv) > 2 else os.path.join(code, "results", "results_default.tsv")
OUT = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(code), "paper", "generated")
os.makedirs(OUT, exist_ok=True)
def out(name): return os.path.join(OUT, name)
# Display order and names. Texts not listed here keep their TSV order, after the listed ones.
ORDER = ["headers", "dna-rand", "dna-rep", "english", "dna", "proteins", "sources", "dblp.xml"]
PIZZA = {"english", "dna", "proteins", "sources", "dblp.xml"}
TEXTNAME = {"dblp.xml": "\\textsf{xml}"}
PER_FLOAT = 3
def order(rows):
    seen = list(collections.OrderedDict((r["text"], 1) for r in rows))
    key = {t: (ORDER.index(t) if t in ORDER else len(ORDER) + i) for i, t in enumerate(seen)}
    return sorted(rows, key=lambda r: key[r["text"]])   # stable: keeps structure order within a text
same, dflt = order(load(a_same)), order(load(a_dflt))
def texts_of(rows): return list(collections.OrderedDict((r["text"], 1) for r in rows))
def tname(t): return TEXTNAME.get(t, "\\textsf{%s}" % t.replace("_", "\\_"))
def fname(t): return t.replace(".", "-")   # file-name-safe text name for p_*.dat
def sname(s):
    s = s.replace("-rrr", "")
    m = {"sada": "cst-sada", "sct3": "cst-sct3", "sct3-t2": "cst-sct3-t2", "cn": "cst-cn", "cnr": "cst-cnr"}
    if s in m: return "\\textsc{%s}" % m[s]
    if s.startswith("fcst:"): return "\\textsc{fcst}$_{%s}$" % s.split(":")[1]
    if s.startswith("ffcst:"): a, b = s.split(":")[1:]; return "\\textsc{ffcst}$_{%s,%s}$" % (a, b)
    return s
def us(x):
    v = float(x) / 1000
    return "%.2f" % v if v < 1 else ("%.1f" % v if v < 10 else "%.0f" % v)
def table(rows, label, caption, fn):
    ts = texts_of(rows)
    chunks = [ts[i:i + PER_FLOAT] for i in range(0, len(ts), PER_FLOAT)]
    o = []
    for ci, chunk in enumerate(chunks):
        lab = label if ci == 0 else "%s-%d" % (label, ci + 1)
        cap = caption if ci == 0 else "%s (continued)." % caption.split(".")[0]
        o += ["\\begin{table*}[tp]\\centering\\small\n\\caption{%s}\\label{%s}\n" % (cap, lab),
              "\\begin{tabular}{@{}lrrrrrrrr@{}}\\toprule\n",
              "Structure & bits/char & \\textsf{Parent} & \\textsf{SDep} & \\textsf{FChild} & \\textsf{NSibling} & \\textsf{Child} & \\textsf{SLink} & \\textsf{LCA}\\\\\n"]
        for t in chunk:
            rs = [r for r in rows if r["text"] == t]
            o.append("\\midrule\n\\multicolumn{9}{@{}l}{%s ($n=%s$)}\\\\\n" % (tname(t), "{:,}".format(int(rs[0]["n"])).replace(",", "\\,")))
            mb = min(float(r["bpc"]) for r in rs)
            for r in rs:
                b = "%.2f" % float(r["bpc"])
                if float(r["bpc"]) == mb: b = "\\textbf{%s}" % b
                o.append("%s & %s & %s\\\\\n" % (sname(r["s"]), b, " & ".join(us(r[k]) for k in ["par","dep","fch","sib","chi","sl","lca"])))
        o.append("\\bottomrule\\end{tabular}\n\\end{table*}\n")
    open(out(fn), "w").write("".join(o))
CSA = "\\texttt{csa\\_wt<\\allowbreak wt\\_huff<\\allowbreak rrr\\_vector<63>>,\\allowbreak 32,\\allowbreak 64>}"
small = [r for r in same if r["text"] not in PIZZA]
pizza = [r for r in same if r["text"] in PIZZA]
table(small, "tab:same", "Space (bits per character, including the CSA) and average time per operation in microseconds, with every structure built on the same CSA (%s). The smallest structure for each text is in bold." % CSA, "tab_same.tex")
table(pizza, "tab:pizza", "Space (bits per character, including the CSA) and average time per operation in microseconds on the Pizza\\&Chili texts, with every structure built on the same CSA (%s). The smallest structure for each text is in bold." % CSA, "tab_pizza.tex")
table(dflt, "tab:default", "Space and time (microseconds) with each library's default configuration. \\textsc{fcst} and \\textsc{ffcst} here use our stand-alone C library (plain wavelet matrix over the BWT, SA and ISA sampled every 32 positions).", "tab_default.tex")
# Construction: one tabular per text group (8 MB texts, Pizza&Chili) inside a single float.
o = ["\\begin{table}[t]\\centering\\small\n\\caption{Construction time (seconds) and peak resident memory (MB) on the same-CSA configuration.}\\label{tab:build}\n"]
for gi, grp in enumerate([g for g in (small, pizza) if g]):
    texts = texts_of(grp)
    if gi: o.append("\\par\\medskip\n")
    o += ["\\begin{tabular}{@{}l%s@{}}\\toprule\n" % ("rr" * len(texts)),
          "& " + " & ".join("\\multicolumn{2}{c}{%s}" % tname(t) for t in texts) + "\\\\\n",
          "Structure" + " & s & MB" * len(texts) + "\\\\\\midrule\n"]
    for s in [r["s"] for r in grp if r["text"] == texts[0]]:
        cells = []
        for t in texts:
            r = [x for x in grp if x["text"] == t and x["s"] == s][0]
            cells += ["%.1f" % float(r["bs"]), "%.0f" % float(r["mb"])]
        o.append("%s & %s\\\\\n" % (sname(s), " & ".join(cells)))
    o.append("\\bottomrule\\end{tabular}\n")
o.append("\\end{table}\n")
open(out("tab_build.tex"), "w").write("".join(o))
for t in texts_of(same):
    for g in ["other", "fcst", "ffcst"]:
        with open(out("p_%s_%s.dat" % (fname(t), g)), "w") as f:
            f.write("name group bpc par sl fch lca\n")
            for r in same:
                s = r["s"].replace("-rrr", "")
                gg = "ffcst" if s.startswith("ffcst") else ("fcst" if s.startswith("fcst") else "other")
                if r["text"] == t and gg == g:
                    f.write("%s %s %s %.3f %.3f %.3f %.3f\n" % (s.replace(":", "-"), g, r["bpc"], float(r["par"])/1e3, float(r["sl"])/1e3, float(r["fch"])/1e3, float(r["lca"])/1e3))
