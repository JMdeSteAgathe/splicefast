"""
Etape 4b : correction FDR, mise en forme et rapports.

FDR : Benjamini-Hochberg PAR ECHANTILLON sur toutes les lignes testees
(exactement la strategie de FRASER : on cherche les outliers d'un individu
contre la cohorte), puis aggregation par gene (p du meilleur site x nombre de
sites testes du gene, puis BH sur les genes).

Sorties (dans <workdir>/results) :
  hits_all.tsv        tous les evenements significatifs, tous echantillons
  <echantillon>.tsv   un fichier par echantillon
  <echantillon>.bed   memes evenements, pour charger dans IGV
  summary.tsv         nombre d'evenements par echantillon et par type
  SpliceFast_viewer.html  visualiseur : ouvre le dossier results, interroge les echantillons
"""
from __future__ import annotations

import html
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .analyze import block_metric, testable
from .classify import Classifier
from .gtf import Annotation
from .stats import bh
from .store import Store
from .util import Progress, Workdir, log, warn, winpath

JUNCTION_METRICS = ("jaccard", "psi5", "psi3")


def _site_blocks(st: Store) -> np.ndarray:
    sb = np.full(st.NS, -1, dtype=np.int32)
    for b in range(st.n_blocks):
        i, j = st.block(b)
        for arr in (st.don_site[i:j], st.acc_site[i:j]):
            s = arr[arr >= 0]
            if len(s):
                sb[s] = b
    return sb


TIER = {"FDR": 0, "gene": 1, "candidate": 2}


def gate_sample(p: np.ndarray, delta: np.ndarray,
                gene_of_col: Optional[np.ndarray], padj_cutoff: float,
                delta_cutoff: float, top_per_sample: int = 20,
                nominal_p: float = 1e-4, top_n: int = 300) -> dict:
    """Selection des lignes d'UN echantillon pour UNE metrique.

    Fonction unique, utilisee a la fois par le rapport et par `why` : le
    diagnostic ne peut plus dire \"rapporte\" pour une ligne que le rapport a
    ecartee (c'etait le cas en 2.3).

    Trois portes, par ordre de priorite :
      FDR      : padj (BH par echantillon) <= seuil et |delta| >= seuil ;
      gene     : a la FRASER2 (`padjustGene`) -- p du gene = min des p des
                 jonctions de TAILLE D'EFFET SUFFISANTE x nombre de jonctions
                 testees du gene, puis BH entre genes. Une seule jonction par
                 gene est ajoutee : la meilleure. En 2.3, le p du gene pouvait
                 venir d'une jonction a delta minuscule (AHNAK : p = 1e-63 pour
                 |delta| < 0.1) et TOUTES les jonctions du gene a |delta| >= 0.1
                 entraient alors, quelle que soit leur p (jusqu'a p = 0.95) ;
      candidat : les `top_per_sample` meilleures p restantes avec p <= nominal_p.

    Plafond `top_n` : applique par ordre de priorite (FDR, puis gene, puis
    candidat) puis de p croissante. En 2.3 le plafond s'appliquait AUX SEULES
    lignes FDR, apres avoir reserve la place des portes secondaires : avec
    > 300 lignes \"gene\", toutes les lignes FDR de la metrique disparaissaient
    (TBCK, SMC1A/jaccard, TUBGCP2/jaccard...).
    """
    ncol = len(p)
    ok = np.isfinite(p)
    padj = bh(p)
    with np.errstate(invalid="ignore"):
        big = ok & (np.abs(delta) >= delta_cutoff)
    reason = np.full(ncol, "", dtype="<U12")   # "candidate" = 9 car. (tronque en <U8)
    with np.errstate(invalid="ignore"):
        reason[big & (padj <= padj_cutoff)] = "FDR"
    gcol = np.full(ncol, np.nan)
    if gene_of_col is not None and ncol:
        g = np.asarray(gene_of_col)
        ng = int(g.max()) + 1 if len(g) else 0
        if ng > 0:
            tested = ok & (g >= 0)
            cnt = np.bincount(g[tested], minlength=ng).astype(np.float64)
            eff = big & (g >= 0)
            gp = np.ones(ng)
            np.minimum.at(gp, g[eff], p[eff])
            gp = np.clip(gp * np.maximum(cnt, 1), 0, 1)
            gp[cnt == 0] = np.nan
            gpadj = bh(gp)
            gcol[tested] = gpadj[g[tested]]
            # meilleure jonction (taille d'effet suffisante) de chaque gene
            ie = np.where(eff)[0]
            if len(ie):
                o = ie[np.lexsort((p[ie], g[ie]))]
                first = np.ones(len(o), dtype=bool)
                first[1:] = g[o][1:] != g[o][:-1]
                best = o[first]
                with np.errstate(invalid="ignore"):
                    best = best[(gcol[best] <= padj_cutoff) & (reason[best] == "")]
                reason[best] = "gene"
    if top_per_sample > 0:
        extra = np.where(big & (reason == "") & (p <= nominal_p))[0]
        if len(extra):
            extra = extra[np.argsort(p[extra])[:top_per_sample]]
            reason[extra] = "candidate"
    cols = np.where(reason != "")[0]
    n_trunc = 0
    if len(cols) > top_n:
        tier = np.array([TIER[r] for r in reason[cols]])
        o = np.lexsort((p[cols], tier))
        n_trunc = len(cols) - top_n
        cols = np.sort(cols[o[:top_n]])
    return dict(cols=cols, reason=reason, padj=padj, gpadj_col=gcol,
                n_trunc=n_trunc, big=big)


def collect_hits(wd: Workdir, st: Store, meta: dict, padj_cutoff: float,
                 delta_cutoff: float, top_n: int,
                 samples: Optional[Sequence[str]] = None,
                 top_per_sample: int = 20,
                 nominal_p: float = 1e-4,
                 noise: Optional[dict] = None) -> List[dict]:
    """Evenements de jonction retenus pour relecture (voir gate_sample).

    `noise` (dict) recoit, par echantillon et par metrique, le nombre de tests
    et le nombre de p <= 1e-4 : c'est la base de l'indice de bruit du resume.
    """
    S = st.S
    sel = (list(range(S)) if not samples
           else [st.samples.index(s) for s in samples if s in st.samples])
    hits: List[dict] = []
    trunc: Dict[str, int] = defaultdict(int)
    for metric in meta["metrics"]:
        ncol = meta["info"][metric]["n_cols"]
        pvp = wd.store / f"pv_{metric}.dat"
        if not pvp.exists():
            continue
        pv = np.memmap(winpath(pvp), dtype=np.float32, mode="r", shape=(S, ncol))
        dl = np.memmap(winpath(wd.store / f"delta_{metric}.dat"),
                       dtype=np.float32, mode="r", shape=(S, ncol))
        gene_of_col = (st.gene if metric != "theta" else None)
        for si in sel:
            row = np.asarray(pv[si], dtype=np.float64)
            ok = np.isfinite(row)
            if not ok.any():
                continue
            p = np.full(ncol, np.nan)
            p[ok] = np.power(10.0, -row[ok])
            delta = np.asarray(dl[si], dtype=np.float64)
            if noise is not None:
                nt = int(ok.sum())
                with np.errstate(invalid="ignore"):
                    ns_ = int(((row >= 4.0) & (np.abs(delta) >= delta_cutoff)).sum())
                noise.setdefault(st.samples[si], {})[metric] = (nt, ns_)
            g = gate_sample(p, delta, gene_of_col, padj_cutoff, delta_cutoff,
                            top_per_sample, nominal_p, top_n)
            if g["n_trunc"]:
                trunc[metric] += g["n_trunc"]
            for c in g["cols"].tolist():
                h = dict(sample=st.samples[si], sample_idx=si, metric=metric,
                         reference=("yes" if st.reference[si] else "no"),
                         retenu_par=str(g["reason"][c]),
                         col=c, p=float(p[c]), padj=float(g["padj"][c]),
                         mlog10p=float(row[c]), delta=float(delta[c]))
                if gene_of_col is not None:
                    h["gene_idx"] = int(gene_of_col[c])
                    h["padj_gene"] = float(g["gpadj_col"][c])
                hits.append(h)
        del pv, dl
    for metric, n in trunc.items():
        warn(f"{metric} : {n} ligne(s) au-dela de --top {top_n} ecartee(s) "
             f"(les moins significatives des portes les moins prioritaires)")
    log(f"{len(hits)} evenements candidats (padj<={padj_cutoff}, "
        f"|delta|>={delta_cutoff})")
    return hits


def cohort_hist(values: np.ndarray, own: float, expected: float,
                own_n: int, bins: int = 25) -> str:
    """Distribution compacte de psi dans la cohorte de reference (JSON).

    {"h": [25 effectifs sur 0..1], "v": psi de l'echantillon, "e": attendu,
     "n": nb d'echantillons de reference couverts, "q": [q05, mediane, q95]}
    Sert au mini-graphique du visualiseur.
    """
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    h, _ = np.histogram(np.clip(v, 0, 1), bins=bins, range=(0.0, 1.0))
    q = ([round(float(x), 4) for x in np.quantile(v, [0.05, 0.5, 0.95])]
         if len(v) else [])
    return json.dumps(dict(h=h.tolist(), v=round(float(own), 4),
                           e=round(float(expected), 4) if np.isfinite(expected) else None,
                           n=int(len(v)), q=q, ns=own_n),
                      separators=(",", ":"))


def enrich_hits(st: Store, hits: List[dict], meta: dict,
                site_block: Optional[np.ndarray], ann: Annotation) -> None:
    """Ajoute k, n, psi, contexte cohorte, en ne lisant chaque bloc qu'une fois."""
    by_metric: Dict[str, Dict[int, List[dict]]] = defaultdict(lambda: defaultdict(list))
    for h in hits:
        c = h["col"]
        if h["metric"] == "theta":
            b = int(site_block[c]) if site_block is not None else -1
        else:
            b = int(st.block_id[c])
        if b >= 0:
            by_metric[h["metric"]][b].append(h)
    for metric, blocks in by_metric.items():
        pr = Progress(len(blocks), f"details {metric}")
        for b, hs in blocks.items():
            pr.update()
            res = block_metric(st, b, metric)
            if res is None:
                continue
            K, N, cols = res
            pos = {int(c): i for i, c in enumerate(cols.tolist())}
            psi = K / np.maximum(N, 1)
            for h in hs:
                i = pos.get(h["col"])
                if i is None:
                    continue
                s = h["sample_idx"]
                h["k"] = int(K[s, i])
                h["n"] = int(N[s, i])
                h["psi"] = float(psi[s, i])
                h["mu"] = float(psi[s, i] - h["delta"])
                col = psi[:, i]
                cov = (N[:, i] >= max(5, meta.get("min_cov", 5))) & st.reference
                h["psi_median_cohorte"] = float(np.median(col[cov])) if cov.any() else float("nan")
                h["n_samples_couverts"] = int(cov.sum())
                others = cov.copy()
                others[s] = False
                h["cohort_hist"] = cohort_hist(col[others], psi[s, i], h["mu"],
                                               int(N[s, i]))
        pr.close()


class ExtensionIndex:
    """Exon allonge / raccourci : jonction dont UN seul site est annote.

    Si le site nouveau tombe DANS l'intron annote qui partage l'autre site,
    l'exon voisin est prolonge dans l'intron : c'est une retention PARTIELLE
    de l'intron (cas NIPBL c.5808+159C>T : donneur cryptique a +157, exon 31
    prolonge de 157 pb). S'il tombe dans l'exon, l'exon est raccourci.
    """

    def __init__(self, ann: Annotation):
        self.by_end: Dict[int, List[int]] = defaultdict(list)
        self.by_start: Dict[int, List[int]] = defaultdict(list)
        for a, b in zip(np.asarray(ann.intron_start).tolist(),
                        np.asarray(ann.intron_end).tolist()):
            self.by_end[int(b)].append(int(a))
            self.by_start[int(a)].append(int(b))

    def describe(self, gs: int, ge: int, start_known: bool, end_known: bool,
                 known_intron: bool, chrom: str = "", local_off: int = 0) -> str:
        if known_intron or start_known == end_known:
            return ""
        if end_known:            # debut nouveau, fin annotee
            cands = self.by_end.get(int(ge), [])
            if not cands:
                return ""
            a = min(cands, key=lambda x: abs(x - gs))
            lo, hi = a, ge
            if a < gs < ge:
                n = gs - a
                return (f"exon extension: +{n} bp into intron {chrom}:"
                        f"{lo - local_off + 1}-{hi - local_off} (partial intron retention)")
            if gs < a:
                return f"exon shortened by {a - gs} bp"
            return ""
        cands = self.by_start.get(int(gs), [])
        if not cands:
            return ""
        b = min(cands, key=lambda x: abs(x - ge))
        if gs < ge < b:
            n = b - ge
            return (f"exon extension: +{n} bp into intron {chrom}:"
                    f"{gs - local_off + 1}-{b - local_off} (partial intron retention)")
        if ge > b:
            return f"exon shortened by {ge - b} bp"
        return ""


def annotate_hits(st: Store, hits: List[dict], ann: Annotation,
                  ref_off: np.ndarray, shift_window: int = 20,
                  pseudo_exon_max: int = 1000,
                  canon: Optional[np.ndarray] = None) -> None:
    cl = Classifier(ann, st.ref_names, ref_off, shift_window=shift_window,
                    pseudo_exon_max=pseudo_exon_max)
    ext = ExtensionIndex(ann)
    # introns / sites des transcrits Ensembl_canonical
    canon_keys = canon if canon is not None else np.zeros(0, np.int64)
    canon_sites = np.unique(np.concatenate([
        canon_keys >> 24, (canon_keys >> 24) + (canon_keys & ((1 << 24) - 1))])) \
        if len(canon_keys) else np.zeros(0, np.int64)

    def _in(arr, v):
        if not len(arr):
            return False
        i = np.searchsorted(arr, v)
        return bool(i < len(arr) and arr[i] == v)

    for h in hits:
        c = h["col"]
        if h["metric"] == "theta":
            gpos = int(ann.site_pos[c]) if ann.n_sites else 0
            kind = int(ann.site_kind[c]) if ann.n_sites else 0
            gi = int(ann.site_gene[c]) if ann.n_sites else -1
            chrom, lp = cl.local(gpos)
            h["chrom"] = chrom
            h["start"] = lp + 1
            h["end"] = lp + 1
            h["gstart"] = gpos
            h["gend"] = None
            h["site_type"] = {1: "donor", 2: "acceptor", 3: "donor+acceptor"}.get(kind, "?")
            h["type"] = ("intron_retention" if h["delta"] < 0
                         else "efficiency_gain")
            if canon is not None:
                h["canonical"] = "yes" if _in(canon_sites, gpos) else "no"
                if h["type"] == "efficiency_gain":
                    h["type"] = ("efficiency_gain_canonical" if h["canonical"] == "yes"
                                 else "efficiency_gain_non-canonical")
            h["gene_idx"] = gi
            h["don_known"] = True
            h["acc_known"] = True
            h["known_intron"] = True
            h["n_skipped_exons"] = 0
            h["comment"] = (
                f"reads a cheval exon/intron en excès au site {h['site_type']}"
                if h["delta"] < 0 else "")
        else:
            gs = int(st.jstart[c])
            ge = int(st.jend[c])
            chrom, ls = cl.local(gs)
            _, le = cl.local(ge)
            h["chrom"] = chrom
            h["start"] = ls + 1        # 1-based, 1re base intronique
            h["end"] = le              # derniere base intronique
            h["gstart"] = gs
            h["gend"] = ge
            h["don_known"] = bool(st.don_site[c] >= 0)
            h["acc_known"] = bool(st.acc_site[c] >= 0)
            h["known_intron"] = bool(st.known_intron[c])
            h["site_type"] = ""
            if canon is not None:
                h["canonical"] = "yes" if _in(canon_keys, int(st.jkey[c])) else "no"
            info = cl.classify_junction(gs, ge, h["don_known"], h["acc_known"],
                                        h["known_intron"], h.get("gene_idx", -1))
            h.update(info)
            h.setdefault("comment", "")
        gname, gid, strand = cl.gene_label(h.get("gene_idx", -1))
        h["gene"] = gname
        h["gene_id"] = gid
        h["strand"] = strand
        if h["metric"] != "theta":
            # 2.5 : le classement donneur/accepteur etait GENOMIQUE (debut /
            # fin de la jonction) ; sur le brin - c'est l'inverse.
            if strand == "-" and h.get("type") in ("novel_donor", "novel_acceptor"):
                h["type"] = ("novel_acceptor" if h["type"] == "novel_donor"
                             else "novel_donor")
            det = ext.describe(h["gstart"], h["gend"], h["don_known"],
                               h["acc_known"], h["known_intron"], chrom=h["chrom"],
                               local_off=h["gstart"] - (h["start"] - 1))
            if det:
                h["event_detail"] = det
        elif strand == "-" and h.get("site_type") in ("donor", "acceptor"):
            h["site_type"] = "acceptor" if h["site_type"] == "donor" else "donor"
        h["intron_len"] = (h["gend"] - h["gstart"]) if h.get("gend") else 0
    # appariement des pseudo-exons, echantillon par echantillon
    per_sample: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        per_sample[h["sample"]].append(h)
    for sid, hs in per_sample.items():
        cl.pair_pseudo_exons(hs)


def merge_metrics(hits: List[dict]) -> List[dict]:
    """Fusionne les lignes identiques detectees par plusieurs metriques."""
    best: Dict[Tuple, dict] = {}
    for h in hits:
        key = (h["sample"], h["chrom"], h["start"], h["end"],
               "theta" if h["metric"] == "theta" else "junction")
        cur = best.get(key)
        if cur is None:
            h["metrics"] = h["metric"]
            best[key] = h
        else:
            cur["metrics"] = ",".join(sorted(set(
                cur["metrics"].split(",") + [h["metric"]])))
            if h["mlog10p"] > cur["mlog10p"]:
                h["metrics"] = cur["metrics"]
                best[key] = h
    out = list(best.values())
    out.sort(key=lambda x: (x["sample"], -x["mlog10p"]))
    return out


# PanelApp (panels, mode de transmission) en tete : c'est ce qu'on lit en
# premier pour trier les anomalies ; les phenotypes, longs, restent a droite
PANEL_COLUMNS = [
    ("panels", "panels"), ("panels_green", "panels_green"),
    ("panel_confidence", "panelapp_best_confidence"), ("moi", "moi"),
]

COLUMNS = PANEL_COLUMNS + [
    ("sample", "sample"), ("run", "run"), ("reference", "in_reference"),
    ("gene", "gene"), ("gene_id", "gene_id"), ("strand", "strand"),
    ("coordinates", "coordinates"), ("chrom", "chrom"), ("start", "start"),
    ("end", "end"),
    ("anomaly_type", "anomaly_type"), ("type", "event_type"),
    ("event_detail", "event_detail"), ("metrics", "metrics"),
    ("k", "split_reads"), ("n", "total_reads"), ("psi", "psi"),
    ("mu", "psi_expected"), ("delta", "delta_psi"),
    ("psi_median_cohorte", "psi_cohort_median"),
    ("p", "p_value"), ("padj", "padj_sample"), ("padj_gene", "padj_gene"),
    ("retenu_par", "selected_by"),
    ("anomalie", "anomaly_id"), ("p_anomalie", "anomaly_p"),
    ("n_preuves", "n_evidence"), ("preuves", "evidence"),
    ("p_combinee", "retention_p_combined"),
    ("padj_combinee", "retention_padj_combined"),
    ("detail_combinaison", "retention_detail"),
    ("n_carriers", "n_carriers"), ("n_carriers_noisy", "n_carriers_noisy"),
    ("carriers", "other_carriers"),
    ("recurrent", "recurrent"), ("hypervariable_locus", "hypervariable_locus"),
    ("intron_len", "intron_length"), ("n_skipped_exons", "skipped_exons"),
    ("don_known", "start_site_annotated"), ("acc_known", "end_site_annotated"),
    ("known_intron", "intron_annotated"), ("canonical", "ensembl_canonical"),
    ("dist_annot", "dist_to_annotated_site"),
    ("site_type", "site_type"), ("n_samples_couverts", "n_ref_samples_covered"),
    ("phenotypes", "phenotypes"),
]
# colonnes JSON : seulement dans les TSV par echantillon (visualiseur)
JSON_COLUMNS = [("panel_details", "panel_details"),
                ("carriers_detail", "carriers_detail"),
                ("cohort_hist", "cohort_hist"), ("context", "context")]


ANOM_COLUMNS = PANEL_COLUMNS + [
    ("sample", "sample"), ("run", "run"), ("gene", "gene"),
    ("coordinates", "coordinates"), ("anomalie", "anomaly_id"),
    ("type", "anomaly_type"), ("p_anomalie", "anomaly_p"),
    ("n_preuves", "n_evidence"), ("preuves", "evidence"),
    ("p_combinee", "retention_p_combined"),
    ("padj_combinee", "retention_padj_combined"),
    ("detail_combinaison", "retention_detail"), ("n_lignes", "n_rows"),
    ("meilleure_ligne", "best_row"), ("retenu_par", "selected_by"),
    ("n_carriers", "n_carriers"), ("recurrent", "recurrent"),
    ("hypervariable_locus", "hypervariable_locus"),
    ("reference", "in_reference"), ("phenotypes", "phenotypes"),
]


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        if math.isnan(v):
            return ""
        if v != 0 and abs(v) < 1e-3:
            return f"{v:.3e}"
        return f"{v:.4g}"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    return str(v)


def _clean(v: str) -> str:
    """Un champ TSV ne doit contenir ni tabulation ni retour a la ligne."""
    return v.replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_tsv(path: Path, hits: List[dict], columns=None) -> None:
    """TSV cite a la maniere d'Excel : les libelles PanelApp contiennent des
    guillemets et des virgules qui, sans citation, decalaient les colonnes."""
    import csv
    # pas de BOM : le classeur .xlsx est la pour Excel, les TSV restent
    # lisibles tels quels par pandas, R ou awk
    with open(winpath(path), "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", quoting=csv.QUOTE_MINIMAL,
                       lineterminator="\n")
        columns = columns or COLUMNS
        w.writerow([h for _, h in columns])
        for h in hits:
            w.writerow([_clean(_fmt(h.get(k))) for k, _ in columns])


def write_xlsx_report(path: Path, hits: List[dict], counts, samples, qc,
                      params: dict, anomalies: Optional[List[dict]] = None) -> None:
    """Classeur Excel : anomalies, evenements, resume, QC, parametres."""
    from .xlsxout import write_xlsx

    def cell(h, k):
        v = h.get(k)
        if isinstance(v, (float, int)) and not isinstance(v, bool):
            return v
        return _clean(_fmt(v))

    types = sorted({h.get("type", "") for h in hits})
    res_rows = [[s, sum(counts.get(s, {}).values())]
                + [counts.get(s, {}).get(t, 0) for t in types] for s in samples]
    qc_rows = [[q.get("sample"), q.get("reads"), q.get("spliced_reads"),
                q.get("junction_reads")] for q in qc]
    sheets = []
    if anomalies:
        sheets.append(("anomalies", [f for _, f in ANOM_COLUMNS],
                       [[cell(a, k) for k, _ in ANOM_COLUMNS] for a in anomalies]))
    sheets += [("events", [f for _, f in COLUMNS],
                [[cell(h, k) for k, _ in COLUMNS] for h in hits]),
               ("summary", ["sample", "total"] + types, res_rows),
               ("qc", ["sample", "reads", "spliced_reads", "junction_reads"], qc_rows),
               ("parameters", ["parameter", "value"],
                [[k, _clean(str(v))] for k, v in params.items()])]
    write_xlsx(path, sheets)


def write_bed(path: Path, hits: List[dict]) -> None:
    with open(winpath(path), "w", encoding="utf-8") as fh:
        fh.write('track name="splicefast" description="evenements aberrants" '
                 'itemRgb="On"\n')
        colors = {"intron_retention": "204,0,0", "exon_skipping": "0,102,204",
                  "pseudo_exon": "153,0,204", "novel_donor": "255,153,0",
                  "novel_acceptor": "255,153,0",
                  "splice_site_shift": "0,153,102"}
        for h in hits:
            s = max(0, int(h["start"]) - 1)
            e = int(h["end"]) if h["end"] else s + 1
            name = f"{h.get('gene','')}|{h.get('type','')}|padj={_fmt(h.get('padj'))}"
            fh.write(f"{h['chrom']}\t{s}\t{max(e, s+1)}\t{name}\t"
                     f"{min(1000, int(10*h['mlog10p']))}\t{h.get('strand','.')}\t"
                     f"{s}\t{max(e, s+1)}\t{colors.get(h.get('type'),'100,100,100')}\n")


def noise_table(noise: Dict[str, dict], cov_noise: Dict[str, Tuple[int, int]]) -> Dict[str, dict]:
    """Indice de bruit par echantillon.

    Taux de tests tres significatifs (p <= 1e-4 et |delta| >= seuil) sur
    tous les tests de jonction/site de l'echantillon, rapporte a la mediane de
    la cohorte. 1 = typique ; > 2 = bruyant ; > 4 = tres bruyant (ARN degrade,
    autre tissu, profondeur faible, contamination...). Les evenements d'un
    echantillon bruyant sont a interpreter avec prudence, et l'echantillon
    devrait sans doute etre sorti de la reference (`reference --test-only`).
    """
    rate = {}
    for sid, d in noise.items():
        nt = sum(v[0] for v in d.values())
        ns = sum(v[1] for v in d.values())
        rate[sid] = (ns + 0.5) / nt if nt else float("nan")
    crate = {sid: ((v[1] + 0.5) / v[0] if v[0] else float("nan")) for sid, v in cov_noise.items()}
    out = {}
    med = np.nanmedian(list(rate.values())) if rate else float("nan")
    cmed = np.nanmedian(list(crate.values())) if crate else float("nan")
    for sid in set(rate) | set(crate):
        r = rate.get(sid, float("nan"))
        c = crate.get(sid, float("nan"))
        ij = r / med if med and np.isfinite(med) and med > 0 else float("nan")
        ic = c / cmed if cmed and np.isfinite(cmed) and cmed > 0 else float("nan")
        ix = np.nanmax([ij, ic]) if np.isfinite([ij, ic]).any() else float("nan")
        level = ("" if not np.isfinite(ix) else "very high" if ix > 4 else
                 "high" if ix > 2 else "low" if ix < 0.5 else "normal")
        out[sid] = dict(noise_index=round(float(ix), 2) if np.isfinite(ix) else None,
                        noise_index_junctions=round(float(ij), 2) if np.isfinite(ij) else None,
                        noise_index_coverage=round(float(ic), 2) if np.isfinite(ic) else None,
                        noise_level=level)
    return out


SUMMARY_COLUMNS = ["sample", "file", "run", "in_reference", "reads", "spliced_reads",
                   "spliced_fraction", "noise_index", "noise_level",
                   "noise_index_junctions", "noise_index_coverage",
                   "n_rows", "n_anomalies", "n_anomalies_fdr",
                   "n_anomalies_panel_green", "n_recurrent_rows"]


def write_summary(path: Path, hits: List[dict], samples: Sequence[str],
                  qc: Sequence[dict], runs: Dict[str, str], ref: np.ndarray,
                  noise: Dict[str, dict]) -> Dict[str, Dict[str, int]]:
    from .util import safe_name
    types = sorted({h.get("type", "?") for h in hits})
    counts: Dict[str, Dict[str, int]] = {s: {t: 0 for t in types} for s in samples}
    anoms: Dict[str, set] = defaultdict(set)
    anoms_fdr: Dict[str, set] = defaultdict(set)
    anoms_green: Dict[str, set] = defaultdict(set)
    nrec: Dict[str, int] = defaultdict(int)
    nrow: Dict[str, int] = defaultdict(int)
    for h in hits:
        sid = h["sample"]
        counts.setdefault(sid, {t: 0 for t in types})
        counts[sid][h.get("type", "?")] = counts[sid].get(h.get("type", "?"), 0) + 1
        nrow[sid] += 1
        a = h.get("anomalie") or f"{h['chrom']}:{h['start']}"
        anoms[sid].add(a)
        if h.get("retenu_par") == "FDR":
            anoms_fdr[sid].add(a)
        if h.get("panels_green"):
            anoms_green[sid].add(a)
        nrec[sid] += h.get("recurrent") == "yes"
    qcm = {q["sample"]: q for q in qc}
    with open(winpath(path), "w", encoding="utf-8") as fh:
        fh.write("\t".join(SUMMARY_COLUMNS + [f"rows_{t}" for t in types]) + "\n")
        for i, s in enumerate(samples):
            q = qcm.get(s, {})
            nz = noise.get(s, {})
            reads = q.get("reads", "")
            sp = q.get("spliced_reads", "")
            frac = (f"{int(sp) / int(reads):.3f}" if reads and sp and int(reads) else "")
            vals = [s, f"{safe_name(s)}.tsv" if nrow.get(s) else "", runs.get(s, ""),
                    "yes" if ref[i] else "no", reads, sp, frac,
                    nz.get("noise_index"), nz.get("noise_level", ""),
                    nz.get("noise_index_junctions"), nz.get("noise_index_coverage"),
                    nrow.get(s, 0), len(anoms.get(s, ())), len(anoms_fdr.get(s, ())),
                    len(anoms_green.get(s, ())), nrec.get(s, 0)]
            vals += [counts.get(s, {}).get(t, 0) for t in types]
            fh.write("\t".join("" if v is None else str(v) for v in vals) + "\n")
    return counts


def add_context(st: Store, hits: List[dict], panelapp=None, legacy=None) -> Dict[str, str]:
    """Run de l'echantillon, coordonnees en une colonne, annotation PanelApp."""
    from .panels import infer_run

    bams = st.meta.get("bams", [])
    runs = {s: infer_run(bams[i] if i < len(bams) else "")
            for i, s in enumerate(st.samples)}
    override = st.wd.root / "runs.tsv" if hasattr(st, "wd") else None
    if override is not None and override.exists():
        n = 0
        for line in open(winpath(str(override)), encoding="utf-8"):
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] in runs:
                runs[f[0]] = f[1].strip()
                n += 1
        log(f"runs.tsv : {n} echantillon(s) avec un run impose")
    cache: Dict[Tuple[str, str], dict] = {}
    n_hit = 0
    for h in hits:
        h["run"] = runs.get(h["sample"], "")
        s_, e_ = h.get("start"), h.get("end")
        h["coordinates"] = (f"{h.get('chrom', '')}:{s_}-{e_}" if e_ and e_ != s_
                            else f"{h.get('chrom', '')}:{s_}")
        key = (str(h.get("gene", "")), str(h.get("gene_id", "")))
        if key not in cache:
            if panelapp is not None:
                cache[key] = panelapp.annotate(*key)
            elif legacy is not None:
                a = legacy.annotate(key[0])
                cache[key] = dict(panels=a.get("panels", ""), panels_green="",
                                  panel_confidence=a.get("panel_conf", ""),
                                  moi=a.get("moi", ""), phenotypes=a.get("phenotypes", ""),
                                  panel_details="")
            else:
                cache[key] = {}
        h.update(cache[key])
        n_hit += bool(h.get("panels"))
    if panelapp is not None or legacy is not None:
        log(f"{n_hit} ligne(s) dans un gene d'au moins un panel")
    return runs


# --------------------------------------------------------------------------- #
# 2.5 : recurrence dans la cohorte, loci hypervariables, contexte exonique
# --------------------------------------------------------------------------- #
# 2.7 : uniquement les genes eux-memes (HLA, immunoglobulines, recepteurs T,
# KIR, LILR, MICA/B). La 2.5-2.6 marquait aussi toute la REGION du CMH
# (chr6:28,5-33,5 Mb) : ~200 genes ordinaires y etaient masques par defaut
# (SYNGAP1, COL11A2, DAXX, BRD2...).
HYPERVARIABLE_GENE = re.compile(
    r"^(HLA-|IG[HKL][VDJC]\d|IGH[AEGMD]\d?$|IG[KL]C\d?$|IGLL|IG[HKL]$|TR[ABDG][VDJC]\d?$|"
    r"TR[ABDG][VDJC]\d|KIR\d|KIR[23]D|"
    r"LILR[AB]\d|MICA$|MICB$)")


def hypervariable(h: dict) -> bool:
    return bool(HYPERVARIABLE_GENE.match(str(h.get("gene", ""))))


def add_recurrence(wd: Workdir, st: Store, meta: dict, hits: List[dict],
                   cov_model: Optional[dict], delta_cutoff: float,
                   p_carrier: float = 1e-4, recurrent_min: int = 4,
                   noisy: Optional[set] = None,
                   site_block: Optional[np.ndarray] = None) -> None:
    """Combien d'echantillons de la cohorte portent le MEME evenement ?

    2.7 - porteur = meme ligne (jonction, site ou intron) avec
      * p <= `p_carrier` (1e-4),
      * ecart dans le MEME sens,
      * |delta| >= max(seuil delta, moitie du delta de l'echantillon) : un
        echantillon a +0.08 n'est pas porteur d'un evenement a +0.60.
    Les echantillons tres bruyants (indice de bruit > 4 : ARN degrade, autre
    tissu...) sont listes a part et ne comptent pas dans `n_carriers` : leur
    "retention" generalisee faisait passer pour recurrents de vrais evenements.
    Le detail de chaque porteur (delta, k/n, p) part dans `carriers_detail`
    (JSON) pour le visualiseur.
    """
    S = st.S
    noisy = noisy or set()
    lim = -np.log10(p_carrier)
    by_metric: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        m = h.get("metric")
        if m in ("coverage", "combined"):
            by_metric["coverage"].append(h)
        elif m in meta.get("metrics", []):
            by_metric[m].append(h)
    for metric, hs in by_metric.items():
        if metric == "coverage":
            if cov_model is None:
                continue
            cs = list(cov_model["samples"])
            idx = {s_: i for i, s_ in enumerate(cs)}
            refm = np.asarray(cov_model["ref"], dtype=bool)
            for h in hs:
                r = h.get("cov_row")
                if r is None:
                    continue
                car = _carrier_idx(h, cov_model["lp"][r].astype(np.float64) / -np.log(10),
                                   cov_model["delta"][r].astype(np.float64),
                                   cov_model["testable"][r], lim, delta_cutoff,
                                   idx.get(h["sample"]))
                _set_carriers(h, car, cs, noisy, refm,
                              lambda i: (float(cov_model["K"][r, i]), float(cov_model["N"][r, i])),
                              cov_model["lp"][r] / -np.log(10), cov_model["delta"][r])
            continue
        ncol = meta["info"][metric]["n_cols"]
        pv = np.memmap(winpath(wd.store / f"pv_{metric}.dat"), dtype=np.float32,
                       mode="r", shape=(S, ncol))
        dl = np.memmap(winpath(wd.store / f"delta_{metric}.dat"), dtype=np.float32,
                       mode="r", shape=(S, ncol))
        cols = np.array(sorted({int(h["col"]) for h in hs}), dtype=np.int64)
        pos = {c: i for i, c in enumerate(cols.tolist())}
        P = np.empty((S, len(cols)), dtype=np.float32)
        D = np.empty((S, len(cols)), dtype=np.float32)
        step = 20000
        for a in range(0, len(cols), step):
            cc = cols[a:a + step]
            P[:, a:a + len(cc)] = pv[:, cc]
            D[:, a:a + len(cc)] = dl[:, cc]
        del pv, dl
        # porteurs, puis k/n des porteurs (bloc relu seulement si necessaire)
        need: Dict[int, List[Tuple[dict, np.ndarray]]] = defaultdict(list)
        for h in hs:
            j = pos[int(h["col"])]
            car = _carrier_idx(h, P[:, j].astype(np.float64), D[:, j].astype(np.float64),
                               np.isfinite(P[:, j]), lim, delta_cutoff, h.get("sample_idx"))
            h["_car"] = (car, j)
            if len([i for i in car if i != h.get("sample_idx")]):
                bk = (int(site_block[h["col"]]) if metric == "theta" and site_block is not None
                      else int(st.block_id[h["col"]]) if metric != "theta" else -1)
                need[bk].append(h)
        kn: Dict[int, Tuple[np.ndarray, np.ndarray]] = {}
        for bk, hl in need.items():
            if bk < 0:
                continue
            res = block_metric(st, bk, metric)
            if res is None:
                continue
            K, N, bcols = res
            bpos = {int(c): i for i, c in enumerate(bcols.tolist())}
            for h in hl:
                i = bpos.get(int(h["col"]))
                if i is not None:
                    kn[id(h)] = (K[:, i].copy(), N[:, i].copy())
        for h in hs:
            car, j = h.pop("_car")
            KN = kn.get(id(h))
            getkn = ((lambda i, KN=KN: (float(KN[0][i]), float(KN[1][i]))) if KN is not None
                     else (lambda i: (float("nan"), float("nan"))))
            _set_carriers(h, car, st.samples, noisy, st.reference, getkn,
                          P[:, j], D[:, j])
    n_rec = 0
    for h in hits:
        h.setdefault("n_carriers", 1)
        h.setdefault("carriers", "")
        h["recurrent"] = "yes" if h["n_carriers"] >= recurrent_min else "no"
        n_rec += h["recurrent"] == "yes"
        h["hypervariable_locus"] = "yes" if hypervariable(h) else "no"
    log(f"recurrence : {n_rec} ligne(s) portee(s) par >= {recurrent_min} "
        f"echantillons fiables (polymorphisme / artefact probable)")


def _carrier_idx(h, mlog10p, delta, ok, lim, delta_cutoff, self_idx) -> List[int]:
    d0 = float(h.get("delta", 0.0) or 0.0)
    sgn = 1.0 if d0 >= 0 else -1.0
    need = max(delta_cutoff, 0.5 * abs(d0))
    with np.errstate(invalid="ignore"):
        car = np.asarray(ok, dtype=bool) & (mlog10p >= lim) & (sgn * delta >= need)
    if self_idx is not None and 0 <= self_idx < len(car):
        car[self_idx] = True
    return np.where(car)[0].tolist()


def _set_carriers(h, car, names, noisy, ref, getkn, mlog10p, delta) -> None:
    self_name = h["sample"]
    others = [i for i in car if str(names[i]) != self_name]
    good = [i for i in others if str(names[i]) not in noisy]
    bad = [i for i in others if str(names[i]) in noisy]
    h["n_carriers"] = 1 + len(good)
    h["n_carriers_noisy"] = len(bad)
    order = sorted(others, key=lambda i: -float(mlog10p[i]))
    txt = [str(names[i]) for i in order if i in good]
    h["carriers"] = ", ".join(txt[:8]) + (f" (+{len(txt) - 8})" if len(txt) > 8 else "") \
        + (f" [+{len(bad)} noisy]" if bad else "")
    det = []
    for i in order[:20]:
        k, n = getkn(i)
        det.append([str(names[i]), round(float(delta[i]), 3),
                    None if k != k else int(round(k)), None if n != n else int(round(n)),
                    float(f"{10 ** -float(mlog10p[i]):.2e}"),
                    1 if str(names[i]) in noisy else 0, 1 if bool(ref[i]) else 0])
    h["carriers_detail"] = json.dumps(det, separators=(",", ":")) if det else ""


def exon_context(ann: Annotation, ref_off: np.ndarray, ref_names: Sequence[str],
                 hits: List[dict], max_exons: int = 16) -> None:
    """Exons du gene autour de l'evenement (JSON), pour le schema du visualiseur."""
    eg = np.asarray(ann.exon_gene)
    es = np.asarray(ann.exon_start)
    ee = np.asarray(ann.exon_end)
    order = np.argsort(eg, kind="stable")
    starts = np.searchsorted(eg[order], np.arange(len(ann.gene_name) + 1))
    cidx = {n: i for i, n in enumerate(ref_names)}
    for h in hits:
        gi = int(h.get("gene_idx", -1) if h.get("gene_idx") is not None else -1)
        c = cidx.get(str(h.get("chrom")))
        if c is None:
            continue
        off = int(ref_off[c])
        a = int(h["start"]) - 1 + off
        b = int(h.get("end") or h["start"]) + off
        span = max(b - a, 1)
        pad = max(400, int(0.6 * span))
        lo, hi = a - pad, b + pad
        ex = []
        if 0 <= gi < len(ann.gene_name):
            ii = order[starts[gi]:starts[gi + 1]]
            s_, e_ = es[ii], ee[ii]
            m = (e_ >= lo) & (s_ <= hi)
            if not m.any():
                # aucun exon dans la fenetre : les deux plus proches
                dist = np.minimum(np.abs(s_ - a), np.abs(e_ - b))
                m = np.zeros(len(ii), bool)
                m[np.argsort(dist)[:2]] = True
            # toujours au moins l'exon le plus proche de CHAQUE cote : un site
            # theta (gain d'efficacite, retention) se lit entre ses deux exons
            left = np.where(e_ <= a + 1)[0]
            right = np.where(s_ >= b - 1)[0]
            if len(left):
                m[left[np.argmax(e_[left])]] = True
            if len(right):
                m[right[np.argmin(s_[right])]] = True
            sel = np.where(m)[0]
            sel = sel[np.argsort(s_[sel])][:max_exons]
            ex = [[int(s_[k]) - off + 1, int(e_[k]) - off] for k in sel]
        h["context"] = json.dumps(dict(x=ex, s=int(h["start"]),
                                       e=int(h.get("end") or h["start"]),
                                       m=str(h.get("metric", "")),
                                       st=h.get("strand", ".")),
                                  separators=(",", ":"))


def _coverage_extras(cov_model: Optional[dict], hits: List[dict],
                     delta_cutoff: float) -> Dict[str, Tuple[int, int]]:
    """Histogramme de cohorte des lignes couverture + bruit couverture."""
    if cov_model is None:
        return {}
    ref = np.asarray(cov_model["ref"], dtype=bool)
    idx = {s: i for i, s in enumerate(cov_model["samples"])}
    for h in hits:
        r = h.get("cov_row")
        if r is None or h.get("metric") not in ("coverage", "combined"):
            continue
        ci = idx.get(h["sample"])
        if ci is None:
            continue
        ok = cov_model["testable"][r] & ref
        ok[ci] = False
        h["cohort_hist"] = cohort_hist(cov_model["psi"][r][ok], float(cov_model["psi"][r, ci]),
                                       float(cov_model["mu"][r, ci]), int(h.get("k") or 0))
    out = {}
    lim = np.log(1e-4)
    for s, ci in idx.items():
        t = cov_model["testable"][:, ci]
        with np.errstate(invalid="ignore"):
            sig = t & (cov_model["lp"][:, ci] <= lim) & (cov_model["delta"][:, ci] >= delta_cutoff)
        out[s] = (int(t.sum()), int(sig.sum()))
    return out


def write_viewer(results: Path) -> Optional[Path]:
    src = Path(__file__).with_name("viewer.html")
    if not src.exists():
        warn("viewer.html absent du paquet : visualiseur non copie")
        return None
    dst = results / "SpliceFast_viewer.html"
    dst.write_bytes(src.read_bytes())
    return dst


def make_reports(wd: Workdir, ann: Annotation, padj_cutoff: float = 0.05,
                 delta_cutoff: float = 0.1, top_n: int = 300,
                 samples: Optional[Sequence[str]] = None,
                 shift_window: int = 20, pseudo_exon_max: int = 1000,
                 panel_files: Optional[Sequence[str]] = None,
                 max_html_rows: int = 15000,
                 cov_min_border: float = 10.0,
                 cov_min_ratio: float = 0.1,
                 cov_min_samples: int = 25,
                 cov_top_per_sample: int = 20,
                 top_per_sample: int = 20,
                 nominal_p: float = 1e-4,
                 harmonise: bool = True,
                 panelapp_cache: Optional[str] = None,
                 recurrent_min: int = 4,
                 write_excel: bool = True) -> List[dict]:
    from .bam import ref_offsets
    from .util import safe_name

    st = Store(wd)
    meta = json.load(open(winpath(wd.store / "analysis.json"), encoding="utf-8"))
    ref_off = ref_offsets(st.ref_lengths)
    noise: Dict[str, dict] = {}
    hits = collect_hits(wd, st, meta, padj_cutoff, delta_cutoff, top_n, samples,
                        top_per_sample=top_per_sample, nominal_p=nominal_p,
                        noise=noise)
    if not hits:
        warn("aucun evenement significatif : essayer --padj 0.1 ou --delta 0.05")
    site_block = _site_blocks(st) if ("theta" in meta["metrics"] and st.NS) else None
    enrich_hits(st, hits, meta, site_block, ann)
    canon = None
    try:
        pj = json.load(open(winpath(wd.store / "project.json"), encoding="utf-8"))
        gtf = pj.get("gtf")
        if gtf and Path(gtf).exists():
            from .gtf import canonical_introns
            canon = canonical_introns(gtf, st.ref_names, ref_off,
                                      cache=wd.store / "canonical_introns.npz")
        elif (wd.store / "canonical_introns.npz").exists():
            canon = np.load(winpath(wd.store / "canonical_introns.npz"))["keys"]
        else:
            warn(f"GTF introuvable ({gtf}) : pas de distinction canonique / "
                 f"non canonique pour efficiency_gain")
    except Exception as e:  # noqa: BLE001
        warn(f"transcrits canoniques non lus ({e.__class__.__name__}: {e})")
    annotate_hits(st, hits, ann, ref_off, shift_window, pseudo_exon_max,
                  canon=canon)
    hits = merge_metrics(hits)
    from .coverage import coverage_hits, coverage_model
    cov_model = coverage_model(wd, min_border=cov_min_border,
                               min_samples=cov_min_samples)
    hits += coverage_hits(wd, ann, padj_cutoff, delta_cutoff,
                          min_border=cov_min_border, min_samples=cov_min_samples,
                          min_ratio=cov_min_ratio, top_n=top_n,
                          top_per_sample=cov_top_per_sample,
                          samples_keep=samples or list(st.samples),
                          model=cov_model)
    anomalies: List[dict] = []
    if harmonise:
        from .harmonize import harmonize
        hits, anomalies = harmonize(wd, st, ann, cov_model, hits, meta,
                                    padj_cutoff=padj_cutoff,
                                    delta_cutoff=delta_cutoff, top_n=top_n,
                                    samples=samples)
    cov_noise = _coverage_extras(cov_model, hits, delta_cutoff)
    noise_pre = noise_table(noise, cov_noise)
    noisy = {s_ for s_, v in noise_pre.items()
             if v.get("noise_index") and v["noise_index"] > 4}
    add_recurrence(wd, st, meta, hits, cov_model, delta_cutoff,
                   recurrent_min=recurrent_min, noisy=noisy,
                   site_block=site_block)
    del cov_model
    gidx = {n: i for i, n in enumerate(ann.gene_name)}
    for h in hits:
        if h.get("gene_idx") is None and h.get("gene") in gidx:
            h["gene_idx"] = gidx[h["gene"]]
    exon_context(ann, ref_off, st.ref_names, hits)
    # anomalies : recurrence / locus hypervariable de leur meilleure ligne
    best: Dict[Tuple[str, str], dict] = {}
    for h in hits:
        k = (h["sample"], h.get("anomalie", ""))
        if k not in best or h["p"] < best[k]["p"]:
            best[k] = h
    for a in anomalies:
        b = best.get((a["sample"], a["anomalie"]))
        if b is not None:
            for f in ("n_carriers", "recurrent", "hypervariable_locus", "gene_id"):
                a[f] = b.get(f)

    # PanelApp (cache telecharge) ; a defaut, anciens exports TSV (--panel)
    from .panelapp import PanelAppDB, find_cache
    pa = None
    pc = find_cache(wd, panelapp_cache)
    if pc is not None:
        pa = PanelAppDB.load(pc)
    legacy = None
    if pa is None and panel_files:
        from .panels import load_panels
        legacy = load_panels(panel_files)
    if pa is None and legacy is None:
        warn("aucune annotation PanelApp : lancer `splicefast panelapp --update`")
    runs = add_context(st, hits, pa, legacy)
    if anomalies:
        add_context(st, anomalies, pa, legacy)
    noise_tab = noise_table(noise, cov_noise)

    wd.results.mkdir(parents=True, exist_ok=True)
    for old in ("result.html", "summary.html", "resultats.xlsx"):
        p = wd.results / old
        if p.exists():
            try:
                p.unlink()
            except OSError:
                pass
    write_tsv(wd.results / "hits_all.tsv", hits)
    if anomalies:
        write_tsv(wd.results / "anomalies.tsv", anomalies, ANOM_COLUMNS)
    per: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        per[h["sample"]].append(h)
    for sid, hs in per.items():
        write_tsv(wd.results / f"{safe_name(sid)}.tsv", hs, COLUMNS + JSON_COLUMNS)
        write_bed(wd.results / f"{safe_name(sid)}.bed", hs)
    qc = st.meta.get("qc", [])
    counts = write_summary(wd.results / "summary.tsv", hits, st.samples, qc,
                           runs, st.reference, noise_tab)
    if pa is not None:
        pa.write_panel_list(wd.results / "panels.tsv")
    with open(winpath(wd.results / "run_info.json"), "w", encoding="utf-8") as fh:
        from . import __version__
        json.dump(dict(version=__version__, padj=padj_cutoff, delta=delta_cutoff,
                       top=top_n, top_per_sample=top_per_sample, nominal_p=nominal_p,
                       recurrent_min=recurrent_min, metrics=meta["metrics"],
                       n_samples=len(st.samples),
                       n_reference=int(st.reference.sum()),
                       panelapp=(pa.source if pa is not None else
                                 ",".join(panel_files or [])),
                       date=__import__("time").strftime("%Y-%m-%d %H:%M")),
                  fh, indent=1)
    vw = write_viewer(wd.results)
    noisy = sorted((v["noise_index"], s) for s, v in noise_tab.items()
                   if v.get("noise_index") and v["noise_index"] > 4)
    if noisy:
        warn(f"{len(noisy)} echantillon(s) tres bruyant(s) (indice > 4), ex. "
             + ", ".join(f"{s} ({x:.1f})" for x, s in noisy[-5:][::-1])
             + " : a sortir de la reference (`reference --test-only`) si autre "
               "tissu ou ARN degrade")
    if write_excel:
        try:
            write_xlsx_report(wd.results / "results.xlsx", hits, counts,
                              list(st.samples), qc,
                              dict(padj=padj_cutoff, delta=delta_cutoff,
                                   top_n=top_n, top_per_sample=top_per_sample,
                                   nominal_p=nominal_p,
                                   metrics=",".join(meta["metrics"]),
                                   samples=len(st.samples),
                                   test_only=int((~st.reference).sum())),
                              anomalies=anomalies)
        except Exception as e:  # noqa: BLE001
            warn(f"classeur Excel non ecrit ({e.__class__.__name__}: {e})")
    log(f"Rapports ecrits dans {wd.results}")
    if vw:
        log(f"Visualiseur : {vw}  (l'ouvrir, puis choisir le dossier results)")
    return hits
