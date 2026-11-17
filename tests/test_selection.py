"""
Test de regression de la selection des jonctions (bug 2.3, cas TBCK).

Configuration calquee sur ECH01 : 400 genes tres exprimes, chacun avec
une jonction a p minuscule mais delta < 0.1 (AHNAK-like) et des jonctions a
|delta| >= 0.1 non significatives ; plus UNE vraie anomalie (TBCK-like,
p = 1e-27, delta = -0.74) dans un autre gene.

    python3 tests/test_selection.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from splicefast.report import gate_sample          # noqa: E402
from splicefast.stats import bh                     # noqa: E402


def selection_23(p, delta, g, padj_cut=0.05, d_cut=0.1, top_n=300,
                 tps=20, nominal=1e-4):
    """Copie fidele de la logique de report.collect_hits en 2.3."""
    ok = np.isfinite(p)
    padj = bh(p)
    big = ok & (np.abs(delta) >= d_cut)
    ng = int(g.max()) + 1
    okg = ok & (g >= 0)
    gp = np.ones(ng)
    np.minimum.at(gp, g[okg], p[okg])
    cnt = np.bincount(g[okg], minlength=ng).astype(float)
    gp = np.clip(gp * np.maximum(cnt, 1), 0, 1)
    gp[cnt == 0] = np.nan
    gpadj = bh(gp)
    strong = big & (padj <= padj_cut)
    reason = np.where(strong, "FDR", "").astype("<U8")
    gcol = np.full(len(p), np.nan)
    gcol[okg] = gpadj[g[okg]]
    by_gene = big & ~strong & (gcol <= padj_cut)
    reason[by_gene] = "gene"
    strong = strong | by_gene
    extra = np.where(big & ~strong & (p <= nominal))[0]
    if len(extra):
        extra = extra[np.argsort(p[extra])[:tps]]
        reason[extra] = "candidat"
        strong[extra] = True
    cols = np.where(strong)[0]
    if len(cols) > top_n:
        keep_extra = cols[reason[cols] != "FDR"]
        fdr = cols[reason[cols] == "FDR"]
        fdr = fdr[np.argsort(p[fdr])[:max(0, top_n - len(keep_extra))]]
        cols = np.concatenate([fdr, keep_extra])
    return set(cols.tolist()), reason


def scenario(seed=0):
    rng = np.random.default_rng(seed)
    p, d, g = [], [], []
    for gene in range(400):
        p.append(10.0 ** -rng.uniform(20, 60)); d.append(rng.uniform(-0.08, 0.08)); g.append(gene)
        for _ in range(3):           # bruit : |delta| >= 0.1, p non significative
            p.append(rng.uniform(0.01, 0.95)); d.append(rng.choice([-1, 1]) * rng.uniform(0.1, 0.5)); g.append(gene)
    for _ in range(200000):          # le reste du transcriptome
        p.append(rng.uniform()); d.append(rng.normal(0, 0.03)); g.append(int(rng.integers(400, 15000)))
    true = len(p)
    p.append(1.3e-27); d.append(-0.74); g.append(15000)   # TBCK-like
    for _ in range(12):
        p.append(rng.uniform()); d.append(rng.normal(0, 0.03)); g.append(15000)
    return np.array(p), np.array(d), np.array(g), true


def main() -> int:
    p, d, g, true = scenario()
    s23, r23 = selection_23(p, d, g)
    out = gate_sample(p, d, g, 0.05, 0.1, 20, 1e-4, 300)
    s24 = set(out["cols"].tolist())
    junk23 = sum(1 for c in s23 if p[c] > 0.05)
    junk24 = sum(1 for c in s24 if p[c] > 0.05)
    print(f"2.3 : {len(s23)} lignes, dont {junk23} a p > 0.05 ; "
          f"vraie anomalie {'RAPPORTEE' if true in s23 else 'PERDUE'}")
    print(f"2.4 : {len(s24)} lignes, dont {junk24} a p > 0.05 ; "
          f"vraie anomalie {'RAPPORTEE (' + out['reason'][true] + ')' if true in s24 else 'PERDUE'}")
    ok = (true in s24) and junk24 == 0 and (true not in s23)
    print("test selection :", "OK" if ok else "ECHEC")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
