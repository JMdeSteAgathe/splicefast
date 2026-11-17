"""
Etape 2 : construction du "store" (reformatage de la matrice en amont).

C'est ici que se trouve la principale difference avec FRASER : au lieu de
manipuler une dizaine de matrices HDF5 de doubles (K, N, psi, delta, mu, rho,
pvals... => plusieurs dizaines de Go en RAM pour 80 echantillons), on stocke :

  * K.split.dat    : matrice int32 (n_echantillons x n_jonctions), memmap
  * K.nonsplit.dat : matrice int32 (n_echantillons x n_sites), memmap
  * junctions.npz  : metadonnees par jonction (position, gene, site, cluster,
                     bloc, annote ou non)

Disposition en lignes = echantillons : lire "toutes les valeurs d'un bloc de
jonctions" revient a lire n_echantillons petits blocs contigus (~100 Ko), donc
tres rapide meme sur un disque externe, et la RAM utilisee est celle d'un seul
bloc (quelques Mo).

Les jonctions sont regroupees en CLUSTERS (composantes connexes partageant un
site d'epissage) puis en BLOCS : un cluster n'est jamais coupe en deux, donc
les psi5/psi3 sont toujours calculables bloc par bloc, ce qui rend l'analyse
"gene par gene" exacte.
"""
from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .bam import LEN_BITS, LEN_MASK, ref_offsets
from .gtf import Annotation
from .util import Progress, Workdir, human, log, warn, winpath, write_json


def _load_sample(path: Path):
    z = np.load(winpath(path), allow_pickle=True)
    return z["jkey"], z["jcount"], z["nonsplit"], z["meta"]


def union_junctions(wd: Workdir, samples: Sequence[Tuple[str, str]],
                    min_total: int = 5, min_samples: int = 1,
                    max_junctions: int = 20_000_000):
    """Union des jonctions de tous les echantillons + statistiques globales.

    Fusion incrementale par numpy (pas de dictionnaire Python geant) :
    memoire ~ 4 tableaux de la taille de l'union.
    """
    keys = np.zeros(0, dtype=np.int64)
    tot = np.zeros(0, dtype=np.int64)
    mx = np.zeros(0, dtype=np.int32)
    ns = np.zeros(0, dtype=np.int16)
    pr = Progress(len(samples), "union des jonctions")
    for sid, _ in samples:
        p = wd.sample_counts(sid)
        if not p.exists():
            warn(f"comptage manquant pour {sid}")
            pr.update()
            continue
        k, c, _, _ = _load_sample(p)
        if len(k) == 0:
            pr.update()
            continue
        newkeys = np.union1d(keys, k)
        if len(newkeys) > max_junctions:
            warn(f"plus de {max_junctions} jonctions : augmenter --min-total")
        idx_old = np.searchsorted(newkeys, keys)
        ntot = np.zeros(len(newkeys), dtype=np.int64)
        nmx = np.zeros(len(newkeys), dtype=np.int32)
        nns = np.zeros(len(newkeys), dtype=np.int16)
        if len(keys):
            ntot[idx_old] = tot
            nmx[idx_old] = mx
            nns[idx_old] = ns
        idx_new = np.searchsorted(newkeys, k)
        ntot[idx_new] += c.astype(np.int64)
        np.maximum.at(nmx, idx_new, c.astype(np.int32))
        nns[idx_new] += 1
        keys, tot, mx, ns = newkeys, ntot, nmx, nns
        pr.update()
    pr.close()
    log(f"Union brute : {len(keys):,} jonctions distinctes")
    keep = (tot >= min_total) & (ns >= min_samples)
    log(f"Apres filtre (total>={min_total}, echantillons>={min_samples}) : "
        f"{int(keep.sum()):,} jonctions")
    return keys[keep], tot[keep], mx[keep], ns[keep]


def _clusters(jstart: np.ndarray, jend: np.ndarray) -> np.ndarray:
    """Composantes connexes des jonctions partageant un site (donneur/accepteur).

    2.5 : calcul EXACT par composantes connexes d'un graphe (scipy). La 2.4
    propageait des etiquettes pendant au plus 200 iterations : sur une longue
    chaine de jonctions (gene de plusieurs centaines d'exons, ou genes voisins
    relies par des jonctions de lecture continue, frequentes avec 300+
    echantillons), la propagation n'avait pas converge et un meme groupe de
    sites pouvait etre coupe entre deux blocs. psi5/psi3/jaccard/theta etaient
    alors calcules sur un groupe INCOMPLET pour les jonctions concernees.
    """
    n = len(jstart)
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    _, dg = np.unique(jstart, return_inverse=True)
    _, ag = np.unique(jend, return_inverse=True)
    nd = int(dg.max()) + 1
    na = int(ag.max()) + 1
    try:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        # graphe biparti jonction <-> site : noeuds 0..n-1 jonctions,
        # n..n+nd-1 donneurs, n+nd.. accepteurs
        rows = np.concatenate([np.arange(n), np.arange(n)])
        cols = np.concatenate([n + dg, n + nd + ag])
        m = n + nd + na
        g = coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)), shape=(m, m))
        _, lab = connected_components(g, directed=False)
        lab = lab[:n].astype(np.int64)
    except ImportError:  # pragma: no cover - repli union-find pur numpy
        parent = np.arange(nd + na, dtype=np.int64)

        def find(x):
            r = x
            while parent[r] != r:
                r = parent[r]
            while parent[x] != r:
                parent[x], x = r, parent[x]
            return r
        for d, a in zip(dg.tolist(), (ag + nd).tolist()):
            rd, ra = find(d), find(a)
            if rd != ra:
                parent[max(rd, ra)] = min(rd, ra)
        lab = np.array([find(d) for d in dg.tolist()], dtype=np.int64)
    # etiquette = plus petit indice de jonction de la composante (stable)
    first = np.full(int(lab.max()) + 1, n, dtype=np.int64)
    np.minimum.at(first, lab, np.arange(n))
    return first[lab]


def check_blocks(jstart: np.ndarray, jend: np.ndarray,
                 block_id: np.ndarray) -> int:
    """Nombre de sites (donneur ou accepteur) partages entre deux blocs.

    Doit valoir 0 : sinon les psi/theta de ces sites sont faux.
    """
    bad = 0
    for pos in (jstart, jend):
        o = np.argsort(pos, kind="stable")
        p, b = pos[o], block_id[o]
        same = p[1:] == p[:-1]
        bad += int((same & (b[1:] != b[:-1])).sum())
    return bad


def _assign_genes(pos_start: np.ndarray, pos_end: np.ndarray,
                  ann: Annotation) -> np.ndarray:
    """Attribue un gene a chaque jonction (site annote, sinon chevauchement)."""
    n = len(pos_start)
    gene = np.full(n, -1, dtype=np.int32)
    if ann.n_sites:
        for pos in (pos_start, pos_end):
            i = np.searchsorted(ann.site_pos, pos)
            ok = (i < ann.n_sites)
            i_ok = np.where(ok, i, 0)
            hit = ok & (ann.site_pos[i_ok] == pos) & (gene < 0)
            gene[hit] = ann.site_gene[i_ok[hit]]
    # jonctions totalement nouvelles : chevauchement avec un gene
    todo = np.where(gene < 0)[0]
    if len(todo) and len(ann.gene_start):
        o = np.argsort(ann.gene_start)
        gs = ann.gene_start[o]
        ge = ann.gene_end[o]
        cmax = np.maximum.accumulate(ge)
        st = pos_start[todo]
        j = np.searchsorted(gs, st, side="right") - 1
        for t, jj, p in zip(todo.tolist(), j.tolist(), st.tolist()):
            k = jj
            steps = 0
            while k >= 0 and cmax[k] > p and steps < 500:
                if ge[k] > p:
                    gene[t] = o[k]
                    break
                k -= 1
                steps += 1
    return gene


# 2.6 : format du store. Un store d'un format anterieur est reconstruit par
# `all` (comptages perdus au remplissage des matrices avant la 2.6).
STORE_FORMAT = 26


def fill_row(keys: np.ndarray, korder: np.ndarray, k: np.ndarray,
             c: np.ndarray) -> np.ndarray:
    """Ligne de comptage d'un echantillon, alignee sur `keys` (non trie).

    Verifie que chaque read d'une jonction retenue est bien place : la
    moindre perte arrete `build` au lieu de produire une matrice fausse.
    """
    J = len(keys)
    row = np.zeros(J, dtype=np.int32)
    if not len(k):
        return row
    sk = keys[korder]
    i = np.searchsorted(sk, k)
    ok = i < J
    i_ok = np.where(ok, i, 0)
    hit = ok & (sk[i_ok] == k)
    row[korder[i_ok[hit]]] = c[hit]
    expect = int(np.asarray(c, dtype=np.int64)[np.isin(k, keys)].sum())
    if int(row.sum(dtype=np.int64)) != expect:
        raise SystemExit("ERREUR interne : comptages mal places dans la matrice "
                         f"({int(row.sum())} reads au lieu de {expect}).")
    return row


def build_store(wd: Workdir, samples: Sequence[Tuple[str, str]], ann: Annotation,
                ref_names: Sequence[str], ref_lengths: Sequence[int],
                min_total: int = 5, min_samples: int = 1,
                block_size: int = 20000) -> dict:
    ref_off = ref_offsets(ref_lengths)
    keys, tot, mx, nsamp = union_junctions(wd, samples, min_total, min_samples)
    if len(keys) == 0:
        raise SystemExit("Aucune jonction retenue : verifier les BAM / les filtres.")

    jstart = (keys >> LEN_BITS).astype(np.int64)
    jlen = (keys & LEN_MASK).astype(np.int64)
    jend = jstart + jlen

    # chromosome + position locale
    chrom = (np.searchsorted(ref_off, jstart, side="right") - 1).astype(np.int32)
    loc_start = jstart - ref_off[chrom]
    loc_end = jend - ref_off[chrom]

    # annotation des sites
    don_site = np.full(len(keys), -1, dtype=np.int64)
    acc_site = np.full(len(keys), -1, dtype=np.int64)
    if ann.n_sites:
        for arr, out in ((jstart, don_site), (jend, acc_site)):
            i = np.searchsorted(ann.site_pos, arr)
            ok = i < ann.n_sites
            i_ok = np.where(ok, i, 0)
            hit = ok & (ann.site_pos[i_ok] == arr)
            out[hit] = i_ok[hit]
    known_intron = np.zeros(len(keys), dtype=bool)
    if len(ann.intron_start):
        ikey = ann.intron_start * np.int64(1 << LEN_BITS) + np.minimum(
            ann.intron_end - ann.intron_start, LEN_MASK)
        ikey.sort()
        i = np.searchsorted(ikey, keys)
        ok = i < len(ikey)
        i_ok = np.where(ok, i, 0)
        known_intron = ok & (ikey[i_ok] == keys)

    gene = _assign_genes(jstart, jend, ann)

    # clusters puis tri : un cluster reste contigu -> blocs exacts
    lab = _clusters(jstart, jend)
    cmin = np.full(len(keys), np.int64(1) << 62, dtype=np.int64)
    np.minimum.at(cmin, lab, jstart)
    order = np.lexsort((jend, jstart, lab, cmin[lab]))

    keys = keys[order]
    jstart, jend, jlen = jstart[order], jend[order], jlen[order]
    chrom, loc_start, loc_end = chrom[order], loc_start[order], loc_end[order]
    don_site, acc_site = don_site[order], acc_site[order]
    known_intron, gene = known_intron[order], gene[order]
    tot, mx, nsamp = tot[order], mx[order], nsamp[order]
    lab = lab[order]
    _, cluster = np.unique(lab, return_inverse=True)
    cluster = cluster.astype(np.int32)

    # blocs : coupures aux frontieres de cluster, ~block_size lignes
    starts = [0]
    last = 0
    change = np.where(np.diff(cluster) != 0)[0] + 1
    for c in change.tolist():
        if c - last >= block_size:
            starts.append(c)
            last = c
    block_starts = np.array(starts + [len(keys)], dtype=np.int64)
    block_id = np.zeros(len(keys), dtype=np.int32)
    for b in range(len(block_starts) - 1):
        block_id[block_starts[b]:block_starts[b + 1]] = b
    log(f"{len(keys):,} jonctions, {int(cluster.max()) + 1:,} clusters, "
        f"{len(block_starts) - 1:,} blocs")
    nbad = check_blocks(jstart, jend, block_id)
    if nbad:
        raise SystemExit(f"ERREUR interne : {nbad} site(s) d'epissage coupe(s) "
                         f"entre deux blocs. Merci de signaler ce message.")

    # ------------------------------------------------------------------ #
    # matrices memmap
    S = len(samples)
    J = len(keys)
    NS = ann.n_sites
    need = S * (J + NS) * 4
    log(f"Matrices : {S} x {J:,} (coupes) + {S} x {NS:,} (non coupes) "
        f"= {human(need)} sur disque")
    ksplit_path = wd.store / "K.split.dat"
    knons_path = wd.store / "K.nonsplit.dat"
    Ks = np.memmap(winpath(ksplit_path), dtype=np.int32, mode="w+", shape=(S, J))
    Kn = (np.memmap(winpath(knons_path), dtype=np.int32, mode="w+", shape=(S, NS))
          if NS else None)

    # 2.6 : les jonctions sont rangees PAR CLUSTER, donc `keys` n'est plus
    # trie. Les versions < 2.6 cherchaient pourtant les comptages par
    # recherche dichotomique dans `keys` : une partie des jonctions de chaque
    # echantillon n'etait pas retrouvee et recevait 0 read (jonctions
    # canoniques a K = 0, d'ou des psi de 1 sur des jonctions minoritaires).
    korder = np.argsort(keys, kind="stable")
    qc = []
    pr = Progress(S, "remplissage des matrices")
    for si, (sid, _) in enumerate(samples):
        p = wd.sample_counts(sid)
        if not p.exists():
            qc.append((sid, 0, 0, 0))
            pr.update()
            continue
        k, c, nsv, meta = _load_sample(p)
        row = fill_row(keys, korder, k, c)
        Ks[si, :] = row
        if Kn is not None:
            if len(nsv) == NS:
                Kn[si, :] = nsv
            elif len(nsv):
                warn(f"{sid} : {len(nsv)} sites au lieu de {NS} (recomptage requis)")
        qc.append((sid, int(meta[0]), int(meta[1]), int(row.sum())))
        pr.update()
    pr.close()
    Ks.flush()
    del Ks
    if Kn is not None:
        Kn.flush()
        del Kn

    np.savez_compressed(
        winpath(wd.store / "junctions.npz"),
        jkey=keys, jstart=jstart, jend=jend, jlen=jlen, chrom=chrom,
        loc_start=loc_start, loc_end=loc_end, don_site=don_site,
        acc_site=acc_site, known_intron=known_intron, gene=gene,
        cluster=cluster, block_id=block_id, block_starts=block_starts,
        total=tot, maxcount=mx, nsamp=nsamp,
    )
    meta = dict(
        n_samples=S, n_junctions=int(J), n_sites=int(NS),
        samples=[s for s, _ in samples], bams=[b for _, b in samples],
        ref_names=list(ref_names), ref_lengths=[int(x) for x in ref_lengths],
        n_blocks=int(len(block_starts) - 1), block_size=block_size,
        store_format=STORE_FORMAT,
        qc=[{"sample": q[0], "reads": q[1], "spliced_reads": q[2],
             "junction_reads": q[3]} for q in qc],
    )
    write_json(wd.store / "store.json", meta)
    return meta


# --------------------------------------------------------------------------- #
# Echantillons de reference / echantillons "hors reference"
#
# Un echantillon hors reference est teste comme les autres, mais ne participe
# pas a la definition de la norme : ni aux facteurs latents, ni a la valeur
# attendue, ni a la dispersion, ni aux medianes de cohorte. Utile pour un
# sequencage trop peu profond, qui sera repasse : son bruit (en particulier la
# fausse retention intronique) n'a plus a etre absorbe par le modele.
# --------------------------------------------------------------------------- #
REFERENCE_FILE = "reference.json"


def read_reference(wd: Workdir, samples: Sequence[str]) -> np.ndarray:
    """Masque booleen (len(samples)) : True = sert de reference."""
    p = wd.store / REFERENCE_FILE
    mask = np.ones(len(samples), dtype=bool)
    if not p.exists():
        return mask
    try:
        d = __import__("json").load(open(winpath(p), encoding="utf-8"))
    except Exception:
        warn("reference.json illisible : tous les echantillons servent de reference")
        return mask
    out = set(d.get("test_only", []))
    for i, sid in enumerate(samples):
        if sid in out:
            mask[i] = False
    return mask


def write_reference(wd: Workdir, test_only: Sequence[str], criteria: dict) -> None:
    write_json(wd.store / REFERENCE_FILE,
               dict(test_only=sorted(set(test_only)), criteria=criteria))


def set_reference(wd: Workdir, min_reads: int = 0,
                  test_only: Sequence[str] = (),
                  reference_only: Optional[Sequence[str]] = None,
                  match: Optional[str] = None) -> Tuple[List[str], List[str]]:
    """Definit qui sert de reference. Renvoie (reference, hors reference).

    min_reads       : sous ce nombre de reads comptes, l'echantillon est teste
                      mais ne sert pas de reference (0 = pas de seuil)
    test_only       : noms a exclure de la reference
    reference_only  : si fourni, SEULS ces noms servent de reference
    """
    meta = __import__("json").load(
        open(winpath(wd.store / "store.json"), encoding="utf-8"))
    samples = list(meta["samples"])
    reads = {q["sample"]: int(q["reads"]) for q in meta.get("qc", [])}
    out: List[str] = []
    reasons: Dict[str, str] = {}
    for sid in samples:
        if reference_only is not None and sid not in set(reference_only):
            out.append(sid)
            reasons[sid] = "hors liste de reference"
        elif sid in set(test_only):
            out.append(sid)
            reasons[sid] = "exclu a la demande"
        elif match and not re.search(match, sid, flags=re.I):
            out.append(sid)
            reasons[sid] = f"nom sans '{match}'"
        elif min_reads and reads.get(sid, 0) < min_reads:
            out.append(sid)
            reasons[sid] = f"{reads.get(sid, 0):,} reads < {min_reads:,}"
    keep = [s for s in samples if s not in set(out)]
    if len(keep) < 20:
        warn(f"seulement {len(keep)} echantillon(s) de reference : la norme "
             f"devient fragile (30 minimum conseilles)")
    write_reference(wd, out, dict(min_reads=int(min_reads), match=match,
                                  test_only=sorted(set(test_only)),
                                  reference_only=(sorted(set(reference_only))
                                                  if reference_only is not None
                                                  else None),
                                  reasons=reasons))
    log(f"reference : {len(keep)} echantillon(s) definissent la norme, "
        f"{len(out)} teste(s) sans y participer")
    for sid in out[:10]:
        log(f"  hors reference : {sid} ({reasons.get(sid, '')})")
    if len(out) > 10:
        log(f"  ... et {len(out) - 10} autre(s)")
    return keep, out


class Store:
    """Acces en lecture au store (memmap, bloc par bloc)."""

    def __init__(self, wd: Workdir):
        self.wd = wd
        self.meta = __import__("json").load(
            open(winpath(wd.store / "store.json"), encoding="utf-8"))
        z = np.load(winpath(wd.store / "junctions.npz"))
        for k in z.files:
            setattr(self, k, z[k])
        self.S = self.meta["n_samples"]
        self.J = self.meta["n_junctions"]
        self.NS = self.meta["n_sites"]
        self.samples = self.meta["samples"]
        self.ref_names = self.meta["ref_names"]
        self.reference = read_reference(wd, self.samples)   # masque (S,)
        self.n_ref = int(self.reference.sum())
        self.ref_lengths = self.meta["ref_lengths"]
        self.Ks = np.memmap(winpath(wd.store / "K.split.dat"), dtype=np.int32,
                            mode="r", shape=(self.S, self.J))
        self.Kn = (np.memmap(winpath(wd.store / "K.nonsplit.dat"),
                             dtype=np.int32, mode="r", shape=(self.S, self.NS))
                   if self.NS else None)

    @property
    def n_blocks(self) -> int:
        return len(self.block_starts) - 1

    def block(self, b: int) -> Tuple[int, int]:
        return int(self.block_starts[b]), int(self.block_starts[b + 1])

    def split_block(self, b: int) -> np.ndarray:
        i, j = self.block(b)
        return np.asarray(self.Ks[:, i:j], dtype=np.float64)

    def split_rows(self, i: int, j: int) -> np.ndarray:
        return np.asarray(self.Ks[:, i:j], dtype=np.float64)
