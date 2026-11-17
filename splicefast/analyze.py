"""
Etape 3 : modele statistique, bloc par bloc.

Pour chaque metrique (jaccard, psi5, psi3, theta) et chaque bloc de jonctions :
    K, N  -> psi = K/N
    mu    -> attendu par regression robuste sur les facteurs latents
    rho   -> dispersion beta-binomiale par jonction
    p     -> test beta-binomial bilateral

On ecrit sur le disque deux matrices float32 par metrique :
    pv_<metrique>.dat     : -log10(p)   (NaN = non teste)
    delta_<metrique>.dat  : psi - mu    (taille de l'effet)
La correction FDR (Benjamini-Hochberg par echantillon, comme FRASER) est faite
au moment du rapport, en lisant UNE ligne (un echantillon) a la fois.

RAM utilisee : un bloc (n_echantillons x block_size) en float64, soit ~13 Mo
pour 80 echantillons et 20 000 jonctions, quel que soit le nombre total de
jonctions du run.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .stats import (bb_logp_screened, choose_subset, design_matrix,
                    latent_factors, logit, refine_rho, rho_moments,
                    robust_expected, trim_mask)
from .store import Store
from .util import Progress, Workdir, human, log, warn, winpath, write_json

METRICS = ("jaccard", "psi5", "psi3", "theta")
LOG10 = np.log(10.0)


def group_sums(K: np.ndarray, keys: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Somme des colonnes de K (S x m) par groupe defini par keys (m).

    Renvoie (N (S x m) = total du groupe de chaque colonne, index de groupe).
    """
    order = np.argsort(keys, kind="stable")
    ks = keys[order]
    starts = np.concatenate(([0], np.nonzero(np.diff(ks))[0] + 1))
    tot = np.add.reduceat(K[:, order], starts, axis=1)   # (S x G)
    gid = np.zeros(len(keys), dtype=np.int64)
    gid[order] = np.repeat(np.arange(len(starts)), np.diff(
        np.concatenate((starts, [len(keys)]))))
    return tot[:, gid], gid


def block_metric(st: Store, b: int, metric: str):
    """Renvoie (K, N, cols) pour la metrique dans le bloc b.

    cols = indices de colonnes dans la matrice de sortie (jonctions ou sites).
    """
    i, j = st.block(b)
    K = st.split_rows(i, j)                     # (S x m) float64
    jstart = st.jstart[i:j]
    jend = st.jend[i:j]
    if metric in ("psi5", "psi3", "jaccard"):
        N5, _ = group_sums(K, jstart)
        N3, _ = group_sums(K, jend)
    if metric == "psi5":
        return K, N5, np.arange(i, j)
    if metric == "psi3":
        return K, N3, np.arange(i, j)
    if metric == "jaccard":
        ns = np.zeros_like(K)
        if st.Kn is not None:
            ds = st.don_site[i:j]
            as_ = st.acc_site[i:j]
            if (ds >= 0).any():
                cols = np.where(ds >= 0)[0]
                ns[:, cols] += np.asarray(st.Kn[:, ds[cols]], dtype=np.float64)
            if (as_ >= 0).any():
                cols = np.where(as_ >= 0)[0]
                ns[:, cols] += np.asarray(st.Kn[:, as_[cols]], dtype=np.float64)
        N = N5 + N3 + ns - K
        return K, np.maximum(N, K), np.arange(i, j)
    if metric == "theta":
        if st.Kn is None:
            return None
        ds = st.don_site[i:j]
        as_ = st.acc_site[i:j]
        sites = np.unique(np.concatenate([ds[ds >= 0], as_[as_ >= 0]]))
        if len(sites) == 0:
            return None
        pos = {s: k for k, s in enumerate(sites.tolist())}
        split = np.zeros((K.shape[0], len(sites)))
        for arr in (ds, as_):
            sel = np.where(arr >= 0)[0]
            if len(sel) == 0:
                continue
            tgt = np.array([pos[s] for s in arr[sel].tolist()], dtype=np.int64)
            np.add.at(split.T, tgt, K[:, sel].T)
        nonsplit = np.asarray(st.Kn[:, sites], dtype=np.float64)
        return split, split + nonsplit, sites
    raise ValueError(metric)


def testable(K: np.ndarray, N: np.ndarray, min_k_one: int = 10,
             min_n_quantile: int = 10, quantile: float = 0.75,
             ref: Optional[np.ndarray] = None) -> np.ndarray:
    """Filtre d'expression a la FRASER.

    `ref` : c'est la cohorte de REFERENCE qui doit etre couverte (denominateur
    N), sinon il n'y a pas de norme a laquelle comparer. Le numerateur K, lui,
    est regarde sur tous les echantillons.
    """
    Nr = N[ref] if (ref is not None and np.any(ref)) else N
    # le NUMERATEUR peut ne venir que d'un seul echantillon, y compris hors
    # reference : c'est exactement le cas d'une jonction nouvelle propre a un
    # patient, qu'il faut tester.
    ok1 = K.max(axis=0) >= min_k_one
    # le DENOMINATEUR, lui, doit etre couvert dans la reference : sans cela il
    # n'y a pas de norme a laquelle comparer.
    ok2 = np.quantile(Nr, quantile, axis=0) >= min_n_quantile
    ok3 = Nr.max(axis=0) > 0
    return ok1 & ok2 & ok3


def collect_subsample(st: Store, metric: str, n_rows: int, min_n: int,
                      seed: int = 0, ref: Optional[np.ndarray] = None):
    """Sous-echantillonne des lignes bien couvertes pour estimer les facteurs."""
    rng = np.random.default_rng(seed)
    Ks: List[np.ndarray] = []
    Ns: List[np.ndarray] = []
    got = 0
    order = rng.permutation(st.n_blocks)
    for b in order:
        res = block_metric(st, int(b), metric)
        if res is None:
            continue
        K, N, _ = res
        keep = testable(K, N, min_k_one=min_n, min_n_quantile=min_n, ref=ref)
        Nr = N[ref] if ref is not None and np.any(ref) else N
        keep &= (Nr >= min_n).mean(axis=0) > 0.8
        if not keep.any():
            continue
        idx = np.where(keep)[0]
        if len(idx) > 2000:
            idx = rng.choice(idx, 2000, replace=False)
        Ks.append(K[:, idx])
        Ns.append(N[:, idx])
        got += len(idx)
        if got >= n_rows:
            break
    if not Ks:
        return None, None
    return np.concatenate(Ks, axis=1), np.concatenate(Ns, axis=1)


def pick_n_factors(eigvals: np.ndarray, S: int, cap: int = 15) -> int:
    w = eigvals[eigvals > 0]
    if len(w) == 0:
        return 2
    q = int((w > w.mean()).sum())
    return int(np.clip(q, 2, min(cap, max(2, S // 6))))


def analyze(wd: Workdir, metrics: Sequence[str] = METRICS,
            n_factors: Optional[int] = None, sub_rows: int = 30000,
            min_k_one: int = 10, min_n_quantile: int = 10,
            min_cov: int = 5, covariates: Optional[np.ndarray] = None,
            robust_iter: int = 3, seed: int = 0,
            n_trim: Optional[int] = None, screen: float = 0.1) -> dict:
    st = Store(wd)
    S, J, NS = st.S, st.J, st.NS
    ref = st.reference
    log(f"Analyse : {S} echantillons, {J:,} jonctions, {NS:,} sites, "
        f"{st.n_blocks} blocs")
    from .store import STORE_FORMAT, check_blocks
    if int(st.meta.get("store_format", 0)) < STORE_FORMAT:
        warn("matrices construites par une version < 2.6 : une partie des "
             "comptages y manque (jonctions a 0 read). Relancer `build` "
             "(sans recompter), puis `analyze` et `report`.")
    nbad = check_blocks(st.jstart, st.jend, st.block_id)
    if nbad:
        warn(f"{nbad} site(s) d'epissage coupe(s) entre deux blocs : matrice "
             f"construite par une version < 2.5 (regroupement inexact). Les "
             f"psi/theta de ces sites sont faux -> relancer `build` (sans "
             f"recompter) puis `analyze`.")
    if not ref.all():
        log(f"Reference : {int(ref.sum())} echantillon(s) definissent la norme ; "
            f"{int((~ref).sum())} teste(s) sans y participer")
        if ref.sum() < 10:
            warn("moins de 10 echantillons de reference : resultats peu fiables")
    info: Dict[str, dict] = {}

    for metric in metrics:
        ncol = NS if metric == "theta" else J
        if metric == "theta" and st.Kn is None:
            warn("pas de comptage non-coupe : metrique theta ignoree")
            continue
        log(f"--- metrique {metric} ---")
        Ksub, Nsub = collect_subsample(st, metric, sub_rows, max(min_cov, 10),
                                       seed, ref=ref)
        if Ksub is None or Ksub.shape[1] < 20:
            warn(f"{metric} : pas assez de lignes couvertes, metrique ignoree")
            continue
        psi = Ksub / np.maximum(Nsub, 1)
        Z_all, eig = latent_factors(psi, n_factors=15, ref=ref)
        from .stats import LAST_FACTOR_INFO
        for ev_, pr_, si_ in LAST_FACTOR_INFO.get("dropped", []):
            warn(f"{metric} : axe ecarte (valeur propre {ev_:.3g}) porte par "
                 f"~{pr_:.1f} echantillon(s), surtout {st.samples[si_]} : "
                 f"echantillon tres atypique, a verifier (QC)")
        q = n_factors or pick_n_factors(eig, int(ref.sum()))
        Z = Z_all[:, :q]
        X = design_matrix(Z, covariates)
        log(f"{metric} : {Ksub.shape[1]} lignes de reference, "
            f"{q} facteurs latents (plan {X.shape[1]} colonnes)")

        pv = np.memmap(winpath(wd.store / f"pv_{metric}.dat"), dtype=np.float32,
                       mode="w+", shape=(S, ncol))
        dl = np.memmap(winpath(wd.store / f"delta_{metric}.dat"), dtype=np.float32,
                       mode="w+", shape=(S, ncol))
        pv[:] = np.nan
        dl[:] = np.nan
        n_tested = 0
        n_exact = 0
        pr = Progress(st.n_blocks, f"test {metric}")
        for b in range(st.n_blocks):
            res = block_metric(st, b, metric)
            pr.update()
            if res is None:
                continue
            K, N, cols = res
            keep = testable(K, N, min_k_one, min_n_quantile, ref=ref)
            if not keep.any():
                continue
            Kk = K[:, keep]
            Nn = N[:, keep]
            cc = cols[keep]
            KT = Kk.T.copy()          # (m x S)
            NT = Nn.T.copy()
            m = KT.shape[0]
            # seuls les echantillons de reference ET informatifs (N >= min_cov)
            # ajustent le modele. En 2.3, les echantillons a N = 0 (gene non
            # exprime) valaient logit(0.5) = 0 : des qu'ils etaient majoritaires
            # (HLA, IG, genes a expression variable), la mediane tombait a 0,
            # Tukey excluait tous les echantillons informatifs et psi attendu
            # valait 0.5 pour tout le monde -- faux delta de +-0.5 d'un cote,
            # dispersion maximale et outliers reels invisibles de l'autre.
            val = (NT >= max(1, min_cov)) & ref[None, :]
            mu = robust_expected(KT, NT, X, n_iter=robust_iter, valid=val)
            msk = trim_mask(KT, NT, mu, n_trim=n_trim, valid=val)
            r0 = rho_moments(KT, NT, mu, mask=msk)
            rho = refine_rho(KT, NT, mu, r0, mask=msk)
            rho_full = np.repeat(rho[:, None], S, axis=1)
            lp, nex = bb_logp_screened(KT, NT, mu, rho_full, screen=screen)
            lp = lp.reshape(m, S)
            n_exact += nex
            low = NT < min_cov
            lp[low] = np.nan
            psi_obs = KT / np.maximum(NT, 1)
            delta = np.where(low, np.nan, psi_obs - mu)
            pv[:, cc] = (-lp / LOG10).T.astype(np.float32)
            dl[:, cc] = delta.T.astype(np.float32)
            n_tested += m
        pr.close()
        pv.flush()
        dl.flush()
        del pv, dl
        info[metric] = dict(n_tested=int(n_tested), n_exact=int(n_exact),
                            n_factors=int(q),
                            n_cols=int(ncol),
                            eigenvalues=[float(x) for x in eig[:15]])
        log(f"{metric} : {n_tested:,} lignes testees")

    meta = dict(metrics=list(info.keys()), info=info, S=S, J=J, NS=NS,
                n_reference=int(ref.sum()),
                test_only=[s for s, r in zip(st.samples, ref) if not r],
                min_cov=min_cov, min_k_one=min_k_one,
                min_n_quantile=min_n_quantile)
    write_json(wd.store / "analysis.json", meta)
    return meta
