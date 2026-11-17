"""
Pourquoi tel evenement n'est-il pas dans le rapport ?

`splicefast why -w DOSSIER --gene DYNC1H1 --sample ECH06-PAX`

Pour chaque jonction du gene (et chaque site pour theta), on affiche ce que le
modele a reellement calcule pour cet echantillon, puis le VERDICT :

  rapporte          : l'evenement est dans les resultats (et par quelle porte)
  sous le seuil     : teste, mais delta ou p-value insuffisants
  non testee        : la jonction n'a jamais atteint le test -- la raison
                      exacte est donnee (pas assez de reads dans la cohorte de
                      reference, couverture de l'echantillon sous --min-cov,
                      jonction absente de la matrice...)

C'est l'outil a lancer quand un evenement connu par ailleurs (FRASER, RT-PCR,
un autre run) n'apparait pas : il dit a quelle etape precise il a disparu.
"""
from __future__ import annotations

from typing import List, Optional

import numpy as np

from .analyze import block_metric, testable
from .gtf import Annotation
from .report import gate_sample
from .store import Store
from .util import Workdir, log, warn, winpath

METRIC_ORDER = ("jaccard", "psi5", "psi3", "theta")


def _gene_indices(ann: Annotation, gene: str) -> List[int]:
    want = gene.strip().upper()
    return [i for i, g in enumerate(ann.gene_name) if str(g).upper() == want]


def _read_rows(wd: Workdir, metric: str, si: int, ncol: int, S: int):
    pvp = wd.store / f"pv_{metric}.dat"
    if not pvp.exists():
        return None, None
    pv = np.memmap(winpath(pvp), dtype=np.float32, mode="r", shape=(S, ncol))
    dl = np.memmap(winpath(wd.store / f"delta_{metric}.dat"), dtype=np.float32,
                   mode="r", shape=(S, ncol))
    return np.asarray(pv[si], dtype=np.float64), np.asarray(dl[si], dtype=np.float64)


def explain(wd: Workdir, ann: Annotation, gene: Optional[str], sample: str,
            chrom: Optional[str] = None, start: Optional[int] = None,
            end: Optional[int] = None, padj_cutoff: float = 0.05,
            delta_cutoff: float = 0.1, nominal_p: float = 1e-4,
            top_per_sample: int = 20, top_n: int = 300) -> List[dict]:
    from .util import read_json

    st = Store(wd)
    meta = read_json(wd.store / "analysis.json")
    if not meta:
        warn("analyse introuvable : lancer 'analyze'")
        return []
    if sample not in st.samples:
        warn(f"echantillon inconnu : {sample}")
        return []
    si = st.samples.index(sample)
    if not st.reference[si]:
        log(f"{sample} est HORS REFERENCE (teste, mais ne definit pas la norme)")

    gi = _gene_indices(ann, gene) if gene else []
    if gene and not gi:
        warn(f"gene absent de l'annotation : {gene}")
    gset = set(gi)

    # colonnes concernees
    sel_j = np.zeros(st.J, dtype=bool)
    if gset:
        sel_j |= np.isin(st.gene, list(gset))
    if chrom is not None:
        ci = st.ref_names.index(chrom) if chrom in st.ref_names else -1
        if ci >= 0:
            off = int(np.cumsum([0] + [int(x) + 1000 for x in st.ref_lengths])[ci])
            a = off + (start or 0) - 5
            b = off + (end or (start or 0)) + 5
            sel_j |= ((st.loc_start >= a) & (st.loc_start <= b)) | \
                     ((st.loc_end >= a) & (st.loc_end <= b))
    cols_j = np.where(sel_j)[0]
    if not len(cols_j):
        warn("aucune jonction de la matrice ne correspond : soit le gene n'a "
             "aucune jonction observee, soit elle a ete filtree a l'etape "
             "'build' (--min-total / --min-samples)")

    out: List[dict] = []
    for metric in meta["metrics"]:
        ncol = meta["info"][metric]["n_cols"]
        row, delta = _read_rows(wd, metric, si, ncol, st.S)
        if row is None:
            continue
        ok = np.isfinite(row)
        p = np.full(ncol, np.nan)
        p[ok] = np.power(10.0, -row[ok])
        # EXACTEMENT la selection du rapport (meme fonction)
        gate = gate_sample(p, delta, st.gene if metric != "theta" else None,
                           padj_cutoff, delta_cutoff, top_per_sample,
                           nominal_p, top_n)
        padj = gate["padj"]
        selected = set(gate["cols"].tolist())
        # rang de l'echantillon parmi ses propres p-values
        order = np.argsort(np.where(np.isfinite(p), p, np.inf))
        rank = np.empty(ncol, dtype=np.int64)
        rank[order] = np.arange(1, ncol + 1)

        if metric == "theta":
            sp = np.asarray(ann.site_pos)
            sg = np.asarray(ann.site_gene)
            cols = np.where(np.isin(sg, list(gset)))[0] if gset else np.zeros(0, int)
            cols = cols[cols < ncol]
        else:
            cols = cols_j

        for c in cols.tolist():
            b = (int(st.block_id[c]) if metric != "theta" else None)
            k = n = float("nan")
            why = ""
            if metric != "theta" and b is not None and b >= 0:
                res = block_metric(st, b, metric)
                if res is not None:
                    K, N, cc = res
                    j = int(np.where(cc == c)[0][0]) if (cc == c).any() else -1
                    if j >= 0:
                        k, n = float(K[si, j]), float(N[si, j])
                        if not np.isfinite(row[c]):
                            keep = testable(K, N, meta.get("min_k_one", 10),
                                            meta.get("min_n_quantile", 10),
                                            ref=st.reference)
                            if not keep[j]:
                                kmax = float(K.max(axis=0)[j])
                                q = float(np.quantile(N[st.reference], 0.75, axis=0)[j])
                                why = (f"jonction non testee : K max cohorte "
                                       f"{kmax:.0f} (>= {meta.get('min_k_one', 10)} "
                                       f"requis), 3e quartile de N dans la "
                                       f"reference {q:.0f} "
                                       f"(>= {meta.get('min_n_quantile', 10)} requis)")
                            elif n < meta.get("min_cov", 5):
                                why = (f"echantillon non teste ici : N = {n:.0f} "
                                       f"< --min-cov {meta.get('min_cov', 5)}")
                            else:
                                why = "non testee (bloc ecarte)"
            d = float(delta[c])
            rep = ""
            if np.isfinite(row[c]):
                if not why:
                    r_ = str(gate["reason"][c])
                    if c in selected:
                        rep = r_
                    elif r_:
                        why = (f"selectionne ({r_}) mais ecarte par --top "
                               f"{top_n} ({gate['n_trunc']} lignes au-dela)")
                    elif abs(d) < delta_cutoff:
                        why = (f"|delta| = {abs(d):.3f} < --delta {delta_cutoff}")
                    elif np.isfinite(gate["gpadj_col"][c]) and gate["gpadj_col"][c] <= padj_cutoff:
                        why = (f"padj = {padj[c]:.3g} ; gene significatif mais une "
                               f"autre jonction du gene est plus probante")
                    else:
                        why = (f"padj = {padj[c]:.3g} > {padj_cutoff} et "
                               f"{rank[c]}e p-value de l'echantillon "
                               f"(> --top-per-sample {top_per_sample})"
                               if p[c] <= nominal_p else
                               f"padj = {padj[c]:.3g} > {padj_cutoff} et "
                               f"p brute {p[c]:.3g} > --nominal-p {nominal_p}")
            out.append(dict(
                metrique=metric, colonne=int(c),
                position=(f"{st.ref_names[int(st.chrom[c])]}:"
                          f"{int(st.loc_start[c]) + 1}-{int(st.loc_end[c])}"
                          if metric != "theta" else f"site {c}"),
                k=k, n=n, psi=(k / n if n else float("nan")),
                psi_attendu=(k / n - d if n else float("nan")),
                delta=d, p=float(p[c]), padj=float(padj[c]),
                rang_echantillon=int(rank[c]) if np.isfinite(p[c]) else -1,
                verdict=("rapporte (" + rep + ")") if rep else (why or "non testee")))
    return out


def print_explanation(rows: List[dict]) -> None:
    if not rows:
        log("rien a expliquer (aucune jonction correspondante)")
        return
    hdr = f"{'metrique':9s} {'position':32s} {'k':>7s} {'n':>7s} {'psi':>6s} " \
          f"{'attendu':>8s} {'delta':>7s} {'p':>10s} {'padj':>10s}  verdict"
    log(hdr)
    rows = sorted(rows, key=lambda r: (r["metrique"], r["position"]))
    for r in rows:
        log(f"{r['metrique']:9s} {r['position']:32s} {r['k']:7.0f} {r['n']:7.0f} "
            f"{r['psi']:6.3f} {r['psi_attendu']:8.3f} {r['delta']:7.3f} "
            f"{r['p']:10.3g} {r['padj']:10.3g}  {r['verdict']}")
