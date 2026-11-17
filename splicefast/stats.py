"""
Coeur statistique - equivalent "leger" du modele de FRASER.

FRASER : psi observe ~ BetaBinomiale(n, mu, rho) ou mu est predit par un
auto-encodeur (denoising autoencoder) entraine sur tous les echantillons, ce
qui corrige les facteurs de confusion (batch, type cellulaire, RIN, sexe...).

Ici on remplace l'auto-encodeur par :
  1. une decomposition en valeurs singulieres (SVD) des logit(psi) sur un
     sous-echantillon de jonctions bien couvertes -> facteurs latents Z (S x q)
     C'est exactement l'espace latent que l'auto-encodeur apprend, mais obtenu
     en une seule passe et en quelques secondes.
  2. par jonction, une regression ponderee robuste de logit(psi) sur [1, Z]
     -> mu attendu pour CHAQUE echantillon (donc controle des covariables)
  3. une dispersion rho par jonction (moments + raffinement par grille)
  4. un test beta-binomial bilateral -> p-value, puis Benjamini-Hochberg
     par echantillon (comme FRASER) et aggregation par gene.

Toutes les fonctions sont vectorisees sur les lignes (jonctions) : un bloc de
20 000 jonctions x 80 echantillons se traite en quelques secondes et quelques
Mo de RAM.
"""
from __future__ import annotations

from typing import Optional, Tuple

import numpy as np

try:  # scipy est fortement recommande mais on garde un repli
    from scipy.special import betainc, betaln, gammaln
    HAVE_SCIPY = True
except Exception:  # pragma: no cover
    HAVE_SCIPY = False
    from math import lgamma

    def gammaln(x):
        x = np.asarray(x, dtype=np.float64)
        out = np.empty_like(x)
        it = np.nditer([x, out], op_flags=[["readonly"], ["writeonly"]])
        for a, b in it:
            b[...] = lgamma(float(a))
        return out

    def betaln(a, b):
        return gammaln(a) + gammaln(b) - gammaln(np.asarray(a) + np.asarray(b))

    def betainc(a, b, x):  # approximation continue si scipy absent
        raise RuntimeError("scipy requis pour les grands n (pip install --user scipy)")

EPS = 1e-8


def logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, EPS, 1 - EPS)
    return np.log(p / (1 - p))


def expit(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))


# --------------------------------------------------------------------------- #
# facteurs latents (remplace l'auto-encodeur)
# --------------------------------------------------------------------------- #

def choose_subset(K: np.ndarray, N: np.ndarray, n_rows: int = 30000,
                  min_n: int = 20, min_frac: float = 0.9,
                  rng: Optional[np.random.Generator] = None) -> np.ndarray:
    """Selectionne des lignes bien couvertes et variables pour la SVD."""
    S = K.shape[0]
    cov_ok = (N >= min_n).mean(axis=0) >= min_frac
    idx = np.where(cov_ok)[0]
    if len(idx) == 0:
        idx = np.argsort(-N.mean(axis=0))[:n_rows]
        return idx
    psi = K[:, idx] / np.maximum(N[:, idx], 1)
    var = psi.std(axis=0)
    # garder les plus variables mais pas les extremes constants
    keep = np.argsort(-var)[: n_rows * 2]
    if len(keep) > n_rows:
        rng = rng or np.random.default_rng(0)
        keep = rng.choice(keep, size=n_rows, replace=False)
    return idx[keep]


LAST_FACTOR_INFO: dict = {}


def latent_factors(psi: np.ndarray, n_factors: int = 10,
                   weights: Optional[np.ndarray] = None,
                   ref: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
    """SVD des logit(psi) centres par ligne. psi : (S x m). Renvoie (Z, valeurs).

    Z : (S x q) facteurs latents des echantillons (orthonormes).

    `ref` (masque S) : les facteurs sont calcules sur ces echantillons
    seulement ; les autres y sont PROJETES. Un echantillon hors reference ne
    peut donc pas tordre les axes du modele, mais recoit quand meme ses
    coordonnees, indispensables pour calculer sa valeur attendue.
    """
    x = logit(psi.T)                       # (m x S)
    S = x.shape[1]
    if ref is None:
        ref = np.ones(S, dtype=bool)
    ref = np.asarray(ref, dtype=bool)
    if ref.sum() < 3:
        ref = np.ones(S, dtype=bool)
    x = x - np.median(x[:, ref], axis=1, keepdims=True)
    x = np.clip(x, -8, 8)
    # 2.4 : winsorisation par ligne a 4 MAD. Un evenement aberrant (le signal
    # que l'on cherche !) ne doit pas peser dans la definition des axes : sans
    # cela, chacun des premiers facteurs pouvait etre porte par UN echantillon
    # aberrant, et la regression "expliquait" ensuite son anomalie.
    mad = 1.4826 * np.median(np.abs(x[:, ref]), axis=1, keepdims=True)
    lim = 4.0 * np.maximum(mad, 0.1)
    x = np.clip(x, -lim, lim)
    if weights is not None:
        x = x * weights[:, None]
    nref = int(ref.sum())
    q = int(min(n_factors, nref - 1, x.shape[0]))
    q = max(q, 1)
    # matrice de covariance S x S (petite) -> eigh, tres rapide et peu de RAM
    C = (x.T @ x) / max(1, x.shape[0])
    w, V = np.linalg.eigh(C[np.ix_(ref, ref)])
    o = np.argsort(-w)
    # 2.4 : un facteur dont l'energie tient dans une poignee d'echantillons
    # n'est pas un facteur de confusion (lot, sexe, RIN...) mais la signature
    # de ces echantillons : on l'ecarte. Nombre effectif d'echantillons d'un
    # axe v (norme 1) = 1 / sum(v^4) ; seuil max(3, 2 % de la reference).
    pr_min = max(3.0, 0.02 * nref)
    keep, dropped = [], []
    for k in o:
        if w[k] <= 0:
            break
        v = V[:, k]
        pr = 1.0 / max(float(np.sum(v ** 4)), 1e-12)
        if pr < pr_min:
            dropped.append((float(w[k]), pr, int(np.argmax(np.abs(v)))))
            continue
        keep.append(k)
        if len(keep) >= q:
            break
    if not keep:
        keep = list(o[:q])
    LAST_FACTOR_INFO.clear()
    sref = np.where(ref)[0]
    LAST_FACTOR_INFO.update(dropped=[(ev, pr, int(sref[i])) for ev, pr, i in dropped])
    keep = np.array(keep)
    w = w[keep]
    q = len(keep)
    Vr = V[:, keep]                        # (Sref x q)
    Z = np.zeros((S, q), dtype=float)
    Z[ref] = Vr
    if (~ref).any():
        # score d'un echantillon hors reference = projection sur les axes de la
        # reference :  z = (x_t . x_ref) V / (m * valeur propre)
        denom = np.maximum(w[:q], 1e-12)
        Z[~ref] = (C[np.ix_(~ref, ref)] @ Vr) / denom[None, :]
        lim = np.abs(Vr).max() * 5 + 1e-9
        Z[~ref] = np.clip(Z[~ref], -lim, lim)   # pas d'extrapolation delirante
    return Z, w


def design_matrix(Z: np.ndarray, extra: Optional[np.ndarray] = None) -> np.ndarray:
    S = Z.shape[0]
    cols = [np.ones((S, 1))]
    if Z.size:
        cols.append(Z)
    if extra is not None and extra.size:
        cols.append(extra)
    return np.concatenate(cols, axis=1)


def wls_rows(Y: np.ndarray, W: np.ndarray, X: np.ndarray,
             ridge: float = 1e-3, leverage: bool = False):
    """Regressions ponderees ligne par ligne, vectorisees.

    Y, W : (m x S). X : (S x p). Renvoie les valeurs ajustees (m x S), et si
    leverage=True le levier h_ss de chaque point (pour la correction
    "leave-one-out" exacte, cf. residus PRESS).
    Cout : un produit matriciel (m x S)(S x p^2) -> tres rapide via BLAS.
    """
    m, S = Y.shape
    p = X.shape[1]
    XX = (X[:, :, None] * X[:, None, :]).reshape(S, p * p)
    A = (W @ XX).reshape(m, p, p)
    A += np.eye(p)[None, :, :] * ridge * np.maximum(W.sum(axis=1), 1)[:, None, None]
    b = (W * Y) @ X
    try:
        beta = np.linalg.solve(A, b[:, :, None])[:, :, 0]
    except np.linalg.LinAlgError:  # pragma: no cover
        beta = np.stack([np.linalg.lstsq(A[i], b[i], rcond=None)[0]
                         for i in range(m)])
    fit = beta @ X.T
    if not leverage:
        return fit
    Ainv = np.linalg.inv(A).reshape(m, p * p)
    h = W * (Ainv @ XX.T)
    return fit, np.clip(h, 0.0, 0.95)


def robust_expected(K: np.ndarray, N: np.ndarray, X: np.ndarray,
                    n_iter: int = 3, huber: float = 3.0, pseudo: float = 0.5,
                    mad_floor: float = 0.05, loo: bool = True,
                    max_leverage: float = 0.5,
                    valid: Optional[np.ndarray] = None) -> np.ndarray:
    """mu attendu par (jonction, echantillon), robuste aux outliers.

    K, N : (m x S) comptages. X : (S x p) plan d'experience [1, Z, covariables].

    Trois protections pour qu'un echantillon aberrant ne "tire" pas son propre
    modele (probleme classique des auto-encodeurs de FRASER) :
      * initialisation sur la MEDIANE de la ligne (depart robuste),
      * re-ponderation de Huber sur plusieurs iterations,
      * correction leave-one-out exacte par le levier (residu PRESS) : la
        valeur attendue pour un echantillon est celle qu'on obtiendrait en
        l'excluant de l'ajustement.
    """
    K = np.asarray(K, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64)
    Y = logit((K + pseudo) / (N + 2 * pseudo))
    W0 = np.sqrt(np.minimum(N, 500.0)) + 1e-3
    # `valid` : echantillons reellement informatifs sur cette ligne. Les autres
    # (couverture trop faible pour que le rapport veuille dire quelque chose)
    # sont exclus de l'ajustement -- sinon leur bruit fixe la norme.
    if valid is not None:
        valid = np.asarray(valid, dtype=bool)
        W0 = W0 * valid
        Yv = np.where(valid, Y, np.nan)
    else:
        Yv = Y
    with np.errstate(invalid="ignore"):
        ref = np.nanmedian(Yv, axis=1, keepdims=True)
    ref = np.where(np.isfinite(ref), ref, np.nanmedian(Y))
    for it in range(max(1, n_iter)):
        r = Y - ref
        rv = np.where(valid, r, np.nan) if valid is not None else r
        with np.errstate(invalid="ignore"):
            s = (1.4826 * np.nanmedian(
                np.abs(rv - np.nanmedian(rv, axis=1, keepdims=True)),
                axis=1, keepdims=True) + mad_floor)
        s = np.where(np.isfinite(s), s, mad_floor)
        # ponderation de Tukey (redescendante) : les points au-dela de
        # huber*MAD sont EXCLUS (poids exactement nul), sinon un outlier garde
        # assez de poids pour tirer sa propre valeur attendue.
        u = np.clip(np.abs(r) / (huber * s), 0.0, 1.0)
        W = W0 * (1.0 - u ** 2) ** 2
        if valid is not None:
            W = np.where(valid, W, 0.0)
        last = it == max(1, n_iter) - 1
        if last and loo:
            fit, h = wls_rows(Y, W, X, leverage=True)
            # correction leave-one-out uniquement si le point ne domine pas son
            # propre ajustement (sinon la correction amplifie le bruit)
            loo_fit = (fit - h * Y) / np.maximum(1.0 - h, 1e-3)
            fit = np.where(h < max_leverage, loo_fit, fit)
        else:
            fit = wls_rows(Y, W, X)
            ref = fit
    return np.clip(expit(fit), 1e-6, 1 - 1e-6)


# --------------------------------------------------------------------------- #
# beta-binomiale
# --------------------------------------------------------------------------- #

def trim_mask(K: np.ndarray, N: np.ndarray, mu: np.ndarray,
              n_trim: Optional[int] = None,
              valid: Optional[np.ndarray] = None) -> np.ndarray:
    """Masque (m x S) excluant les points les plus aberrants de chaque ligne.

    Indispensable pour ne pas laisser l'outlier gonfler sa propre dispersion :
    avec peu d'echantillons, un seul evenement aberrant suffit sinon a rendre
    la jonction "sur-dispersee" donc non significative.
    """
    S = K.shape[1]
    ok = N > 0
    if valid is not None:   # le rognage ne porte que sur les echantillons
        ok = ok & valid     # qui servent de reference
        S = int(np.max(valid.sum(axis=1))) if valid.ndim == 2 else int(valid.sum())
    if n_trim is None:
        n_trim = max(1, int(np.ceil(0.02 * S)))
    n_trim = int(min(n_trim, max(0, S - 4)))
    if n_trim <= 0:
        return ok
    v = np.maximum(mu * (1 - mu), 1e-6)
    chi = np.where(ok, (K - N * mu) ** 2 / (np.maximum(N, 1) * v), -np.inf)
    thr = -np.sort(-chi, axis=1)[:, n_trim - 1][:, None]
    return ok & (chi < thr)


def rho_moments(K: np.ndarray, N: np.ndarray, mu: np.ndarray,
                min_rho: float = 1e-4, max_rho: float = 0.5,
                mask: Optional[np.ndarray] = None) -> np.ndarray:
    """Dispersion par la methode des moments (Pearson), corrigee des ddl."""
    ok = (N > 0) if mask is None else mask
    v = np.maximum(mu * (1 - mu), 1e-6)
    chi = np.where(ok, (K - N * mu) ** 2 / (np.maximum(N, 1) * v), 0.0)
    n_ok = ok.sum(axis=1)
    chi_sum = chi.sum(axis=1)
    denom = np.where(ok, N - 1.0, 0.0).sum(axis=1)
    rho = (chi_sum - np.maximum(n_ok - 1, 1)) / np.maximum(denom, 1e-9)
    return np.clip(rho, min_rho, max_rho)


def bb_loglik(K, N, mu, rho):
    a = mu * (1 - rho) / rho
    b = (1 - mu) * (1 - rho) / rho
    return (betaln(K + a, N - K + b) - betaln(a, b))


def refine_rho(K, N, mu, rho0, n_grid: int = 7, span: float = 8.0,
               mask: Optional[np.ndarray] = None):
    """Raffine rho par une petite grille log autour de l'estimateur des moments.

    La vraisemblance est calculee en excluant les points aberrants (mask), pour
    la meme raison que dans trim_mask.
    """
    w = np.ones_like(K, dtype=np.float64) if mask is None else mask.astype(float)

    def ll_of(r):
        return (bb_loglik(K, N, mu, r[:, None]) * w).sum(axis=1)

    best = rho0.copy()
    bestll = ll_of(rho0)
    grid = np.exp(np.linspace(-np.log(span), np.log(span), n_grid))
    for g in grid:
        if abs(g - 1.0) < 1e-12:
            continue
        r = np.clip(rho0 * g, 1e-5, 0.6)
        ll = ll_of(r)
        better = ll > bestll
        best[better] = r[better]
        bestll[better] = ll[better]
    return best


def _bb_logpmf(x, n, a, b):
    return (gammaln(n + 1) - gammaln(x + 1) - gammaln(n - x + 1)
            + betaln(x + a, n - x + b) - betaln(a, b))


def _logsumexp(a: np.ndarray, axis: int = -1) -> np.ndarray:
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    return (m + np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True)))[
        (slice(None),) * (a.ndim - 1) + (0,)] if axis in (-1, a.ndim - 1) else None


def bb_logp(K: np.ndarray, N: np.ndarray, mu: np.ndarray, rho: np.ndarray,
            window: int = 512, row_chunk: int = 2000) -> np.ndarray:
    """log de la p-value beta-binomiale bilaterale.

    Strategie (validee contre scipy.stats.betabinom) :
      * somme EXACTE des pmf en espace log, depuis k et en s'eloignant du mode
        (la queue est dominee par le terme en k, donc la somme converge vite) ;
      * si la fenetre n'a pas converge (queue lourde), on prend le maximum
        entre cette somme, l'approximation Beta  P(X<=k) ~ I_q(a,b)  et une
        approximation normale avec la variance beta-binomiale : c'est
        conservateur (jamais de p-value trop petite).
    Travailler en log evite toute perte de precision pour p < 1e-300.
    """
    K = np.asarray(K, dtype=np.float64).ravel()
    N = np.asarray(N, dtype=np.float64).ravel()
    mu = np.asarray(mu, dtype=np.float64).ravel()
    rho = np.asarray(rho, dtype=np.float64).ravel()
    n = len(K)
    a = np.maximum(mu * (1 - rho) / rho, 1e-8)
    b = np.maximum((1 - mu) * (1 - rho) / rho, 1e-8)
    logp = np.zeros(n, dtype=np.float64)
    left = K <= N * mu
    valid = N > 0

    W = max(8, int(window) // 2)
    step = np.arange(W, dtype=np.float64)[None, :]
    idx = np.where(valid)[0]
    for s0 in range(0, len(idx), row_chunk):
        ii = idx[s0:s0 + row_chunk]
        k = np.round(K[ii])[:, None]
        nn = N[ii][:, None]
        aa = a[ii][:, None]
        bb = b[ii][:, None]
        lf = left[ii][:, None]
        lo = np.where(lf, 0.0, k)
        hi = np.where(lf, k, nn)

        def lpmf(x):
            return _bb_logpmf(x, nn, aa, bb)

        # fenetre A : depuis l'extremite la plus proche du mode (hi..hi-W+1)
        xa = hi - step
        oka = (xa >= lo) & (xa >= 0)
        lpa = np.where(oka, lpmf(np.where(oka, xa, 0.0)), -np.inf)
        # fenetre B : a l'extremite opposee (lo..lo+W-1), sans recouvrir A
        xb = lo + step
        okb = (xb <= hi) & (xb <= hi - W) & (xb <= nn)
        lpb = np.where(okb, lpmf(np.where(okb, xb, 0.0)), -np.inf)
        # borne majorante rigoureuse pour l'intervalle non couvert
        glo = lo + W
        ghi = hi - W
        L = ghi - glo + 1.0
        hasgap = L > 0
        gl = lpmf(np.clip(glo, 0, nn))
        gh = lpmf(np.clip(ghi, 0, nn))
        lgap = np.where(hasgap, np.log(np.maximum(L, 1.0)) + np.maximum(gl, gh),
                        -np.inf)

        allp = np.concatenate([lpa, lpb, lgap], axis=1)
        mx = np.max(allp, axis=1, keepdims=True)
        mx = np.where(np.isfinite(mx), mx, 0.0)
        tail = (mx + np.log(np.sum(np.exp(allp - mx), axis=1, keepdims=True)))[:, 0]
        logp[ii] = np.minimum(tail, 0.0)

    logp[~valid] = 0.0
    return np.minimum(logp + np.log(2.0), 0.0)


def _erfc(x):
    try:
        from scipy.special import erfc
        return erfc(x)
    except Exception:  # pragma: no cover
        from math import erfc as _e
        return np.vectorize(_e)(x)


def bb_pvalue(K, N, mu, rho, **kw) -> np.ndarray:
    return np.exp(bb_logp(K, N, mu, rho, **kw))


# --------------------------------------------------------------------------- #
# corrections de tests multiples
# --------------------------------------------------------------------------- #

def bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg ; les NaN restent NaN."""
    p = np.asarray(p, dtype=np.float64)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out
    idx = np.where(ok)[0]
    pv = p[idx]
    o = np.argsort(pv)
    ranked = pv[o] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[o] = np.clip(ranked, 0, 1)
    out[idx] = res
    return out


def zscores(psi: np.ndarray, axis: int = 0) -> np.ndarray:
    """Z-scores robustes sur logit(psi) (comme FRASER)."""
    x = logit(psi)
    med = np.median(x, axis=axis, keepdims=True)
    sd = 1.4826 * np.median(np.abs(x - med), axis=axis, keepdims=True)
    sd = np.maximum(sd, 1e-3)
    return (x - med) / sd


def bb_logp_approx(K, N, mu, rho):
    """Approximation normale (avec correction de continuite) du test BB.

    Sert de PRE-FILTRE : elle est bonne dans le corps de la distribution, ou
    l'on n'a de toute facon pas besoin d'une precision extreme, et permet de
    ne lancer le calcul exact (couteux) que sur les quelques points qui
    peuvent devenir significatifs.
    """
    K = np.asarray(K, dtype=np.float64)
    N = np.asarray(N, dtype=np.float64)
    mean = N * mu
    var = N * mu * (1.0 - mu) * (1.0 + (N - 1.0) * rho)
    sd = np.sqrt(np.maximum(var, 1e-12))
    z = (np.abs(K - mean) - 0.5) / sd
    z = np.maximum(z, 0.0)
    # log(2 * Phi(-z)) calcule de facon stable
    return np.log(2.0) + _log_ndtr_neg(z)


def _log_ndtr_neg(z):
    """log(Phi(-z)) pour z >= 0, stable jusqu'a des z tres grands."""
    out = np.empty_like(z)
    small = z < 8.0
    if small.any():
        with np.errstate(divide="ignore"):
            out[small] = np.log(0.5 * _erfc(z[small] / np.sqrt(2.0)))
    big = ~small
    if big.any():
        zb = z[big]
        # developpement asymptotique de Mills
        out[big] = (-0.5 * zb ** 2 - np.log(zb) - 0.5 * np.log(2 * np.pi)
                    + np.log1p(-1.0 / zb ** 2 + 3.0 / zb ** 4))
    return out


def bb_logp_screened(K, N, mu, rho, screen: float = 0.1, **kw):
    """log p beta-binomial, exact la ou ca compte, approche ailleurs.

    Renvoie (logp, n_exact). Les points dont l'approximation normale donne
    deja p > `screen` ne peuvent pas survivre a la correction FDR : on garde
    la valeur approchee. Gain typique : x15 a x30 sur le temps de calcul.
    """
    K = np.asarray(K, dtype=np.float64).ravel()
    N = np.asarray(N, dtype=np.float64).ravel()
    mu = np.asarray(mu, dtype=np.float64).ravel()
    rho = np.asarray(rho, dtype=np.float64).ravel()
    lp = bb_logp_approx(K, N, mu, rho)
    hit = lp <= np.log(screen)
    n = int(hit.sum())
    if n:
        lp[hit] = bb_logp(K[hit], N[hit], mu[hit], rho[hit], **kw)
    return lp, n
