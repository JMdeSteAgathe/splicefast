"""
Harmonisation : reunir les signaux disparates d'une meme anomalie.

Une retention d'intron laisse jusqu'a QUATRE traces independantes ou presque :

  couverture       le corps de l'intron est couvert         (module couverture, delta > 0)
  theta_donneur    reads non coupes a cheval sur le donneur  (theta, delta < 0)
  theta_accepteur  idem a l'accepteur                        (theta, delta < 0)
  jaccard          la jonction canonique perd sa part        (jaccard, delta < 0)

psi5 et psi3 n'en font PAS partie : ils ne comparent que des reads coupes
entre eux. Un transcrit qui retient l'intron ne produit aucun read coupe ; si
aucune autre jonction ne partage le donneur ou l'accepteur, psi5 = psi3 = 1 chez
le patient comme dans la cohorte. Ils sont aveugles a une retention pure par
construction (c'est le cas SMC1A).

Pour chaque echantillon et chaque intron annote mesure par le module de
couverture, on combine ces quatre p-values UNILATERALES (dans le sens d'une
retention) par la methode de Brown : Fisher corrige de la dependance entre
metriques. La dependance (les memes reads non coupes entrent dans theta et dans
jaccard) est mesuree empiriquement sur la cohorte elle-meme : moyenne, variance
et covariances de -2 ln p sur un sous-echantillon de (echantillon, intron).
Sans cette correction, Fisher compterait plusieurs fois la meme information et
fabriquerait des p-values beaucoup trop petites.

Consequences :
  * trois ou quatre signaux concordants, meme moderes, donnent une p combinee
    tres faible : l'anomalie remonte en tete (SMC1A) ;
  * un signal isole contredit par les autres (couverture seule, theta et
    jaccard normaux : transcription antisens, gene imbrique...) est pondere a
    la baisse : c'est le gain de specificite sur les retentions ;
  * une retention dont aucun signal n'atteint seul le seuil peut etre
    rapportee grace a la concordance (gain de sensibilite).

Les autres evenements (jonctions nouvelles, sauts d'exon...) sont regroupes en
anomalies par sites partages, sans combinaison statistique : leur p d'anomalie
est la meilleure p du groupe.
"""
from __future__ import annotations

from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .bam import LEN_BITS, LEN_MASK, ref_offsets
from .stats import bh
from .util import log, warn, winpath

FAMILIES = ("coverage", "theta_donor", "theta_acceptor", "jaccard")
SIGN = np.array([+1.0, -1.0, -1.0, -1.0])   # sens d'une retention
CAP = 2.0 * np.log(1e6)      # -2 ln p plafonne pour estimer les moments
VAR_INFLATION = 1.10
NOMINAL = 1e-3               # une \"preuve\" = p unilaterale <= 1e-3
ABERRANT_TYPES = {"exon_skipping", "pseudo_exon", "novel_donor",
                  "novel_acceptor", "splice_site_shift", "novel_junction"}


def _log_chi2_sf(y: np.ndarray, f: np.ndarray) -> np.ndarray:
    """log P(chi2_f > y), stable dans les queues extremes."""
    from scipy.special import gammaln
    from scipy.stats import chi2
    y = np.asarray(y, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    with np.errstate(all="ignore"):
        out = chi2.logsf(y, f)
    bad = ~np.isfinite(out) & (y > 0)
    if bad.any():
        a = f[bad] / 2.0
        x = y[bad] / 2.0
        # developpement asymptotique de Gamma(a, x) / Gamma(a) pour x >> a
        out[bad] = ((a - 1.0) * np.log(x) - x - gammaln(a)
                    + np.log1p(np.maximum((a - 1.0) / x, -0.99)))
    return np.minimum(out, 0.0)


def one_sided(p2: np.ndarray, delta: np.ndarray, sign: float) -> np.ndarray:
    """p bilaterale + sens de l'ecart -> p unilaterale dans le sens `sign`."""
    with np.errstate(invalid="ignore"):
        good = sign * delta > 0
        p1 = np.where(good, p2 / 2.0, 1.0 - p2 / 2.0)
    p1[~(np.isfinite(p2) & np.isfinite(delta))] = np.nan
    return np.clip(p1, 1e-320, 1.0)


class RetentionEvidence:
    """Les quatre familles de preuves, alignees sur les introns de couverture."""

    def __init__(self, wd, st, ann, model: dict, meta: dict):
        self.wd, self.st, self.model, self.meta = wd, st, model, meta
        d = model["d"]
        cols = model["cols"]
        self.M = len(cols)
        chrom = np.asarray(d["chrom"])[cols]
        start0 = np.asarray(d["start"], dtype=np.int64)[cols]
        end0 = np.asarray(d["end"], dtype=np.int64)[cols]
        gene = np.asarray(d["gene"], dtype=np.int64)[cols]
        self.chrom = chrom
        self.start1 = start0 + 1          # 1re base intronique (1-based)
        self.end1 = end0                   # derniere base intronique
        self.gene = gene
        ref_off = ref_offsets(st.ref_lengths)
        cidx = {n: i for i, n in enumerate(st.ref_names)}
        off = np.array([ref_off[cidx[c]] if c in cidx else -1 for c in chrom],
                       dtype=np.int64)
        gs = np.where(off >= 0, off + start0, -1)
        ge = np.where(off >= 0, off + end0, -1)
        strand = np.asarray(ann.gene_strand)[np.clip(gene, 0, len(ann.gene_strand) - 1)]
        self.minus = strand < 0
        # colonnes theta : site en debut genomique / en fin genomique
        sp = np.asarray(ann.site_pos)
        self.site_lo = self._lookup(sp, gs)
        self.site_hi = self._lookup(sp, ge)
        # colonne jaccard de la jonction canonique
        key = gs * np.int64(1 << LEN_BITS) + np.minimum(ge - gs, LEN_MASK)
        o = np.argsort(st.jkey)
        jk = np.asarray(st.jkey)[o]
        i = np.searchsorted(jk, key)
        i_ok = np.clip(i, 0, len(jk) - 1)
        hit = (gs >= 0) & (i < len(jk)) & (jk[i_ok] == key)
        self.jcol = np.where(hit, o[i_ok], -1)
        # memmaps
        self.mm = {}
        S = st.S
        for metric in ("theta", "jaccard"):
            if metric in meta.get("metrics", []):
                n = meta["info"][metric]["n_cols"]
                pv = wd.store / f"pv_{metric}.dat"
                dl = wd.store / f"delta_{metric}.dat"
                if pv.exists() and dl.exists():
                    self.mm[metric] = (
                        np.memmap(winpath(pv), dtype=np.float32, mode="r", shape=(S, n)),
                        np.memmap(winpath(dl), dtype=np.float32, mode="r", shape=(S, n)))
        self.cov_index = {s: i for i, s in enumerate(model["samples"])}
        self.c = None
        self.f = None

    @staticmethod
    def _lookup(sorted_pos: np.ndarray, pos: np.ndarray) -> np.ndarray:
        if len(sorted_pos) == 0:
            return np.full(len(pos), -1, dtype=np.int64)
        i = np.searchsorted(sorted_pos, pos)
        i_ok = np.clip(i, 0, len(sorted_pos) - 1)
        hit = (pos >= 0) & (i < len(sorted_pos)) & (sorted_pos[i_ok] == pos)
        return np.where(hit, i_ok, -1)

    def _take(self, metric: str, si: int, cols: np.ndarray):
        out_p = np.full(len(cols), np.nan)
        out_d = np.full(len(cols), np.nan)
        if metric not in self.mm:
            return out_p, out_d
        pv, dl = self.mm[metric]
        ok = cols >= 0
        if ok.any():
            rp = np.asarray(pv[si], dtype=np.float64)
            rd = np.asarray(dl[si], dtype=np.float64)
            out_p[ok] = np.power(10.0, -rp[cols[ok]])
            out_d[ok] = rd[cols[ok]]
        return out_p, out_d

    def sample(self, sid: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        """(P unilaterales M x 4, delta M x 4) pour un echantillon."""
        if sid not in self.st.samples or sid not in self.cov_index:
            return None
        si = self.st.samples.index(sid)
        ci = self.cov_index[sid]
        m = self.model
        P2 = np.full((self.M, 4), np.nan)
        D = np.full((self.M, 4), np.nan)
        lp = m["lp"][:, ci].astype(np.float64)
        p = np.exp(lp)
        p[~m["testable"][:, ci]] = np.nan
        P2[:, 0] = p
        D[:, 0] = m["delta"][:, ci]
        tlo, dlo = self._take("theta", si, self.site_lo)
        thi, dhi = self._take("theta", si, self.site_hi)
        # donneur = debut genomique sur le brin +, fin genomique sur le brin -
        P2[:, 1] = np.where(self.minus, thi, tlo)
        D[:, 1] = np.where(self.minus, dhi, dlo)
        P2[:, 2] = np.where(self.minus, tlo, thi)
        D[:, 2] = np.where(self.minus, dlo, dhi)
        P2[:, 3], D[:, 3] = self._take("jaccard", si, self.jcol)
        P1 = np.column_stack([one_sided(P2[:, k], D[:, k], SIGN[k])
                              for k in range(4)])
        return P1, D

    def fit_moments(self, samples: Sequence[str], n_samples: int = 60,
                    seed: int = 0) -> None:
        """Moments empiriques de -2 ln p (dependance entre familles)."""
        rng = np.random.default_rng(seed)
        pool = [s for s in samples if s in self.cov_index and s in self.st.samples]
        if len(pool) > n_samples:
            pool = list(rng.choice(pool, n_samples, replace=False))
        Xs = []
        for sid in pool:
            r = self.sample(sid)
            if r is None:
                continue
            X = np.minimum(-2.0 * np.log(r[0]), CAP)
            if len(X) > 50000:
                X = X[rng.choice(len(X), 50000, replace=False)]
            Xs.append(X)
        mean = np.full(4, 2.0)
        cov = np.diag(np.full(4, 4.0))
        if Xs:
            X = np.concatenate(Xs)
            for a in range(4):
                xa = X[:, a]
                fa = np.isfinite(xa)
                if fa.sum() >= 200:
                    mean[a] = xa[fa].mean()
                    cov[a, a] = max(xa[fa].var(), 1e-3)
            for a in range(4):
                for b in range(a + 1, 4):
                    both = np.isfinite(X[:, a]) & np.isfinite(X[:, b])
                    if both.sum() >= 200:
                        xa, xb = X[both, a], X[both, b]
                        cab = float(np.mean((xa - xa.mean()) * (xb - xb.mean())))
                        # une covariance negative rendrait le test anti-
                        # conservateur : on l'ignore
                        cov[a, b] = cov[b, a] = max(cab, 0.0)
        self.mean, self.cov = mean, cov
        # table par motif de disponibilite (bitmask 1..15)
        self.c = np.ones(16)
        self.f = np.full(16, 2.0)
        for mask in range(1, 16):
            idx = [k for k in range(4) if mask >> k & 1]
            E = mean[idx].sum()
            V = cov[np.ix_(idx, idx)].sum()
            if len(idx) > 1:
                # l'approximation de Brown est un peu optimiste dans la queue
                # quand les familles sont tres correlees (x1.25 a p = 1e-6 en
                # simulation) ; +10 % de variance la rend conservatrice
                # (x0.55 a x0.75 de 1e-3 a 1e-6, 20 M de tirages nuls)
                V *= VAR_INFLATION
            self.c[mask] = V / (2.0 * E)
            self.f[mask] = 2.0 * E * E / V
        corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
        log("harmonisation : correlations de -2 ln p entre familles "
            + ", ".join(f"{FAMILIES[a]}~{FAMILIES[b]} {corr[a, b]:.2f}"
                        for a in range(4) for b in range(a + 1, 4)))

    def combine(self, P1: np.ndarray):
        """-> (log p combinee, nombre de familles, nombre de preuves nominales)."""
        avail = np.isfinite(P1)
        mask = (avail * (1 << np.arange(4))[None, :]).sum(axis=1)
        X = np.where(avail, -2.0 * np.log(np.where(avail, P1, 1.0)), 0.0).sum(axis=1)
        lp = np.zeros(len(P1))
        ok = mask > 0
        lp[ok] = _log_chi2_sf(X[ok] / self.c[mask[ok]], self.f[mask[ok]])
        n_av = avail.sum(axis=1)
        n_nom = (np.where(avail, P1, 1.0) <= NOMINAL).sum(axis=1)
        return lp, n_av, n_nom


# --------------------------------------------------------------------------- #
def _fmt_p(p: float) -> str:
    return "NA" if not np.isfinite(p) else (f"{p:.1e}" if p < 1e-3 else f"{p:.3f}")


def harmonize(wd, st, ann, model: Optional[dict], hits: List[dict], meta: dict,
              padj_cutoff: float = 0.05, delta_cutoff: float = 0.1,
              top_n: int = 300,
              samples: Optional[Sequence[str]] = None) -> Tuple[List[dict], List[dict]]:
    """Ajoute aux hits : anomalie, p_anomalie, n_preuves, preuves.

    Renvoie (hits completes et re-triees, table des anomalies).
    """
    per_sample: Dict[str, List[dict]] = defaultdict(list)
    for h in hits:
        per_sample[h["sample"]].append(h)

    ev = None
    if model is not None and len(model["cols"]):
        ev = RetentionEvidence(wd, st, ann, model, meta)
        ev.fit_moments(list(st.samples))
    else:
        log("harmonisation : pas de module couverture, regroupement simple")

    # index des introns de couverture
    by_intron: Dict[Tuple[str, int, int], int] = {}
    by_site: Dict[Tuple[str, int], List[int]] = defaultdict(list)
    if ev is not None:
        for r in range(ev.M):
            c = str(ev.chrom[r])
            s, e = int(ev.start1[r]), int(ev.end1[r])
            by_intron[(c, s, e)] = r
            by_site[(c, s)].append(r)
            by_site[(c, e + 1)].append(r)

    def rows_of(h) -> List[int]:
        if h.get("metrics") == "coverage" and "cov_row" in h:
            return [int(h["cov_row"])]
        if h.get("metric") == "theta":
            return list(by_site.get((h["chrom"], int(h["start"])), []))
        r = by_intron.get((h.get("chrom"), int(h.get("start") or 0),
                           int(h.get("end") or 0)))
        return [r] if r is not None else []

    todo = (list(samples) if samples else list(st.samples))
    anomalies: List[dict] = []
    out_hits: List[dict] = []
    n_added = 0
    for sid in todo:
        hs = per_sample.get(sid, [])
        rec: Dict[int, dict] = {}
        if ev is not None:
            r = ev.sample(sid)
            if r is not None:
                P1, D = r
                lp, n_av, n_nom = ev.combine(P1)
                pc = np.exp(lp)
                pc_t = np.where(n_av >= 2, pc, np.nan)
                padj = bh(pc_t)
                with np.errstate(invalid="ignore"):
                    effect = ((D[:, 0] >= delta_cutoff) | (D[:, 1] <= -delta_cutoff)
                              | (D[:, 2] <= -delta_cutoff) | (D[:, 3] <= -delta_cutoff))
                    sig = (padj <= padj_cutoff) & (n_nom >= 2) & effect
                want = set(np.where(sig)[0].tolist())
                if len(want) > top_n:
                    w = np.array(sorted(want))
                    want = set(w[np.argsort(lp[w])[:top_n]].tolist())
                for h in hs:
                    want.update(rows_of(h))
                for rr in want:
                    rec[rr] = dict(row=rr, logp=float(lp[rr]), p=float(pc[rr]),
                                   padj=float(padj[rr]), n_av=int(n_av[rr]),
                                   n_nom=int(n_nom[rr]), sig=bool(sig[rr]),
                                   P1=P1[rr].copy(), D=D[rr].copy())
        # ---- regroupement des hits de l'echantillon (union-find) ----------
        parent = list(range(len(hs)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            a, b = find(i), find(j)
            if a != b:
                parent[b] = a

        site_owner: Dict[Tuple[str, int], int] = {}
        row_owner: Dict[int, int] = {}
        hit_row: Dict[int, Optional[int]] = {}
        for i, h in enumerate(hs):
            rows = [x for x in rows_of(h) if x in rec]
            # un site theta partage par plusieurs introns : le plus probant
            best = min(rows, key=lambda x: rec[x]["logp"]) if rows else None
            hit_row[i] = best
            if best is not None:
                if best in row_owner:
                    union(row_owner[best], i)
                else:
                    row_owner[best] = i
            if h.get("metric") == "theta" or h.get("metrics") == "coverage":
                keys = [(h["chrom"], int(h["start"]))]
                if h.get("metrics") == "coverage":
                    keys = [(h["chrom"], int(h["start"])), (h["chrom"], int(h["end"]) + 1)]
            else:
                keys = [(h["chrom"], int(h["start"])), (h["chrom"], int(h["end"]) + 1)]
            for k in keys:
                if k in site_owner:
                    union(site_owner[k], i)
                else:
                    site_owner[k] = i
        groups: Dict[int, List[int]] = defaultdict(list)
        for i in range(len(hs)):
            groups[find(i)].append(i)

        # retentions significatives sans aucun hit : ligne de synthese
        covered_rows = {r for r in hit_row.values() if r is not None}
        synth = []
        for rr, x in rec.items():
            if x["sig"] and rr not in covered_rows:
                synth.append(rr)
        groups_list = [[hs[i] for i in g] for g in groups.values()]
        rows_list = [[hit_row[i] for i in g if hit_row[i] is not None]
                     for g in groups.values()]
        for rr in synth:
            h = _synthetic_hit(sid, rr, rec[rr], ev, ann, model, st)
            groups_list.append([h])
            rows_list.append([rr])
            n_added += 1

        for gi, (members, rows) in enumerate(zip(groups_list, rows_list)):
            rows = sorted(set(rows), key=lambda x: rec[x]["logp"])
            # unites de preuve : chaque intron combine REMPLACE ses composants
            units = []
            fams = set()
            for rr in rows:
                units.append(rec[rr]["p"])
                for k in range(4):
                    if np.isfinite(rec[rr]["P1"][k]) and rec[rr]["P1"][k] <= NOMINAL:
                        fams.add(FAMILIES[k])
            folded = set()
            for h in members:
                comp = _is_component(h, rows, ev)
                if comp:
                    folded.add(id(h))
                    continue
                units.append(h["p"])
                if h["p"] <= NOMINAL:
                    m = h.get("metrics") or h.get("metric")
                    for mm in str(m).split(","):
                        fams.add(mm if mm != "theta" else "theta")
            p_anom = float(min(units)) if units else float("nan")
            best = min(members, key=lambda h: h["p"])
            gene = next((h.get("gene") for h in members if h.get("gene")), "")
            chrom = best["chrom"]
            span = [h for h in members if h.get("metric") != "theta"] or members
            lo = min(int(h["start"]) for h in span)
            hi = max(int(h.get("end") or h["start"]) for h in span)
            aid = f"{gene or chrom}:{lo}-{hi}"
            ret = rows[0] if rows else None
            detail = ""
            typ = best.get("type", "")
            p_comb = padj_comb = float("nan")
            aberrant = [h for h in members if h.get("type") in ABERRANT_TYPES
                        and h["p"] <= NOMINAL]
            if aberrant:
                typ = min(aberrant, key=lambda h: h["p"])["type"]
            if ret is not None:
                x = rec[ret]
                p_comb, padj_comb = x["p"], x["padj"]
                detail = " | ".join(
                    f"{FAMILIES[k]} p={_fmt_p(x['P1'][k])} d={x['D'][k]:+.2f}"
                    for k in range(4) if np.isfinite(x["P1"][k]))
                if x["sig"] or x["n_nom"] >= 2:
                    # une jonction nouvelle significative dans le meme groupe :
                    # c'est elle qui explique la "retention" (site cryptique
                    # intra-intronique = retention PARTIELLE)
                    typ = (f"{typ}+partial_retention" if aberrant
                           else "intron_retention")
            if any("partial intron retention" in str(h.get("event_detail", ""))
                   for h in members):
                typ = "exon_extension (partial intron retention)"
            for h in members:
                h["anomaly_type"] = typ
                h["anomalie"] = aid
                h["p_anomalie"] = p_anom
                h["n_preuves"] = len(fams)
                h["preuves"] = ",".join(sorted(fams))
                if ret is not None:
                    h["p_combinee"] = p_comb
                    h["padj_combinee"] = padj_comb
                    h["detail_combinaison"] = detail
                out_hits.append(h)
            anomalies.append(dict(
                sample=sid, run="", gene=gene, chrom=chrom, start=lo, end=hi,
                anomalie=aid, type=typ, n_lignes=len(members),
                n_preuves=len(fams), preuves=",".join(sorted(fams)),
                p_anomalie=p_anom, p_combinee=p_comb, padj_combinee=padj_comb,
                detail_combinaison=detail,
                meilleure_ligne=f"{best.get('metrics', best.get('metric'))} "
                                f"{best['chrom']}:{best['start']}-{best.get('end') or best['start']} "
                                f"p={_fmt_p(best['p'])} delta={best.get('delta', float('nan')):+.2f}",
                retenu_par=",".join(sorted({str(h.get("retenu_par", "")) for h in members})),
                reference=best.get("reference", "")))
    # hits d'echantillons non traites (ne devrait pas arriver)
    done = {id(h) for h in out_hits}
    for h in hits:
        if id(h) not in done:
            h.setdefault("p_anomalie", h["p"])
            out_hits.append(h)
    out_hits.sort(key=lambda h: (h["sample"], h.get("p_anomalie", h["p"]),
                                 h.get("anomalie", ""), h["p"]))
    anomalies.sort(key=lambda a: (a["sample"], a["p_anomalie"]))
    if n_added:
        log(f"harmonisation : {n_added} retention(s) rapportee(s) grace a la "
            f"concordance des metriques (aucune ne passait seule)")
    log(f"harmonisation : {len(anomalies)} anomalie(s) pour {len(out_hits)} ligne(s)")
    return out_hits, anomalies


def _is_component(h: dict, rows: List[int], ev) -> bool:
    """Le hit est-il deja compte dans la p combinee d'un intron du groupe ?"""
    if ev is None or not rows:
        return False
    m = h.get("metrics") or h.get("metric")
    d = h.get("delta", float("nan"))
    if not np.isfinite(d):
        return False
    # un ecart de sens oppose (gain d'efficacite) n'est PAS une preuve de
    # retention : il reste une unite de preuve a part entiere
    if (m == "coverage" and d <= 0) or (m != "coverage" and d >= 0):
        return False
    for rr in rows:
        c = str(ev.chrom[rr])
        s, e = int(ev.start1[rr]), int(ev.end1[rr])
        if h.get("chrom") != c:
            continue
        if m == "coverage" and int(h["start"]) == s and int(h["end"]) == e:
            return True
        if h.get("metric") == "theta" and int(h["start"]) in (s, e + 1):
            return True
        if (h.get("metric") == "jaccard" and int(h["start"]) == s
                and int(h.get("end") or 0) == e):
            return True
    return False


def _synthetic_hit(sid, rr, x, ev, ann, model, st) -> dict:
    ci = ev.cov_index[sid]
    g = int(ev.gene[rr])
    gname = ann.gene_name[g] if 0 <= g < len(ann.gene_name) else ""
    gid = ann.gene_id[g] if 0 <= g < len(ann.gene_id) else ""
    strand = "-" if ev.minus[rr] else "+"
    return dict(
        sample=sid, gene=gname, gene_id=gid, chrom=str(ev.chrom[rr]),
        start=int(ev.start1[rr]), end=int(ev.end1[rr]), strand=strand,
        type="intron_retention", metrics="combined", metric="combined",
        k=int(round(float(model["K"][rr, ci]))), n=int(round(float(model["N"][rr, ci]))),
        psi=round(float(model["psi"][rr, ci]), 4),
        mu=round(float(model["mu"][rr, ci]), 4),
        delta=round(float(model["delta"][rr, ci]), 4),
        psi_median_cohorte=round(float(model["med"][rr]), 4),
        p=x["p"], padj=x["padj"], mlog10p=float(-x["logp"] / np.log(10)),
        padj_gene=float("nan"), retenu_par="concordance",
        intron_len=int(ev.end1[rr] - ev.start1[rr] + 1),
        reference=("yes" if model["ref"][ci] else "no"),
        n_samples_couverts=int(model["n_cov"][rr]), cov_row=rr)
