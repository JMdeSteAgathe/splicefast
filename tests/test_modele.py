"""
Tests de regression du modele (2.4).

1. psi attendu = 0.5 : gene exprime chez 40 % des echantillons seulement.
2. facteurs latents porteurs d'un seul echantillon aberrant.

    python3 tests/test_modele.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from splicefast.stats import (bb_logp, design_matrix, latent_factors,  # noqa: E402
                              refine_rho, rho_moments, robust_expected,
                              trim_mask)


def test_mu_half() -> bool:
    rng = np.random.default_rng(1)
    S = 100
    N = np.where(np.arange(S) < 60, 0, rng.integers(60, 200, S)).astype(float)
    K = rng.binomial(N.astype(int), 0.03).astype(float)
    K[99] = 0.6 * N[99]
    KT, NT = K[None, :], N[None, :]
    valid = NT >= 5                       # ce que fait analyze en 2.4
    mu = robust_expected(KT, NT, np.ones((S, 1)), valid=valid)
    msk = trim_mask(KT, NT, mu, valid=valid)
    rho = refine_rho(KT, NT, mu, rho_moments(KT, NT, mu, mask=msk), mask=msk)
    p = float(np.exp(bb_logp(KT[:, 99:], NT[:, 99:], mu[:, 99:], rho)[0]))
    ok = abs(mu[0, 70] - 0.03) < 0.02 and p < 1e-20
    print(f"[1] psi attendu {mu[0, 70]:.3f} (vrai 0.03), outlier p = {p:.1e} -> "
          f"{'OK' if ok else 'ECHEC'}")
    return ok


def test_factors() -> bool:
    rng = np.random.default_rng(2)
    S, m = 60, 3000
    batch = np.where(np.arange(S) % 2 == 0, 0.15, -0.15)
    base = rng.normal(1.5, 1.0, m)
    load = rng.normal(0, 1, m)
    x = base[None, :] + batch[:, None] * load[None, :] + rng.normal(0, 0.15, (S, m))
    x[7, :1500] -= 3.0                     # echantillon 7 : aberrant sur la moitie
                                           # du transcriptome (ARN degrade...)
    psi = 1 / (1 + np.exp(-x))
    Z, _ = latent_factors(psi, n_factors=3)
    X = design_matrix(Z)
    # une jonction ou l'echantillon 7 est aberrant
    N = np.full((1, S), 200.0)
    K = np.round(N * psi[:, 0][None, :])
    mu = robust_expected(K, N, X)
    gap = abs(mu[0, 7] - np.median(psi[:, 0]))
    corr = abs(np.corrcoef(Z[:, 0], batch)[0, 1])
    ok = gap < 0.15 and corr > 0.8
    print(f"[2] 1er facteur ~ lot : |r| = {corr:.2f} ; attendu ech. 7 = {mu[0, 7]:.3f} "
          f"(cohorte {np.median(psi[:, 0]):.3f}, observe {psi[7, 0]:.3f}) -> "
          f"{'OK' if ok else 'ECHEC'}")
    return ok


def test_clusters() -> bool:
    """2.5 : regroupement exact meme sur une tres longue chaine de jonctions."""
    from splicefast.store import _clusters
    ok = True
    for nex in (400, 800, 3000):
        s, e = [], []
        for i in range(nex - 1):
            s.append(i * 100 + 10); e.append((i + 1) * 100)
            if i + 2 < nex:
                s.append(i * 100 + 10); e.append((i + 2) * 100)
        o = np.random.default_rng(0).permutation(len(s))
        n = len(np.unique(_clusters(np.array(s)[o], np.array(e)[o])))
        ok &= n == 1
    print(f"[3] chaine de 3000 exons : 1 seul groupe -> {'OK' if ok else 'ECHEC'}")
    return ok


def test_fill() -> bool:
    """2.6 : comptages places dans la matrice malgre l'ordre PAR CLUSTER.

    Gene X : jonctions 1000->9000 et 5000->9000 (accepteur partage, cluster
    de debut 1000) ; gene Y imbrique : 3000->3500 (cluster de debut 3000).
    Ordre par cluster : X(1000), X(5000), Y(3000) -> cles non triees.
    """
    from splicefast.store import fill_row
    js = np.array([1000, 5000, 3000]); je = np.array([9000, 9000, 3500])
    keys = (js << 24) | (je - js)
    k = np.sort(keys); c = np.array([1142, 10, 50])[np.argsort(keys)]
    row = fill_row(keys, np.argsort(keys), k, c)
    i = np.searchsorted(keys, k)                  # ancien code (< 2.6)
    ok = i < len(keys); io = np.where(ok, i, 0); hit = ok & (keys[io] == k)
    old = np.zeros(3, int); old[io[hit]] = c[hit]
    good = row.sum() == c.sum() and old.sum() < c.sum()
    print(f"[4] remplissage matrice : 2.6 {row.tolist()} / avant {old.tolist()} -> "
          f"{'OK' if good else 'ECHEC'}")
    return bool(good)


if __name__ == "__main__":
    r = [test_mu_half(), test_factors(), test_clusters(), test_fill()]
    sys.exit(0 if all(r) else 1)
