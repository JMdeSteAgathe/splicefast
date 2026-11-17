"""
Lecture de l'annotation GTF (GENCODE / Ensembl / RefSeq) et construction :
  * la liste triee des sites d'epissage annotes (coordonnees globales)
  * les introns annotes
  * les exons (pour detecter les saut d'exon)
  * les genes (pour l'aggregation et le rapport)

Tout est stocke en tableaux numpy et mis en cache dans un .npz : le GTF n'est
lu qu'une seule fois par projet.
"""
from __future__ import annotations

from array import array
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .util import log, open_maybe_gzip, warn, winpath

ATTR_KEYS = ("gene_id", "gene_name", "transcript_id", "gene_type", "gene_biotype")


def _attr(s: str, key: str) -> Optional[str]:
    i = s.find(key)
    while i != -1:
        # doit etre precede d'un debut de champ
        if i == 0 or s[i - 1] in " ;\t":
            j = i + len(key)
            while j < len(s) and s[j] in ' =':
                j += 1
            if j < len(s) and s[j] == '"':
                k = s.find('"', j + 1)
                return s[j + 1:k]
            k = j
            while k < len(s) and s[k] not in ';\n':
                k += 1
            return s[j:k].strip().strip('"')
        i = s.find(key, i + 1)
    return None


def build_ref_map(gtf_chroms: Sequence[str], ref_names: Sequence[str]) -> Dict[str, int]:
    """Fait correspondre les noms de chromosomes du GTF a ceux du BAM."""
    idx = {n: i for i, n in enumerate(ref_names)}
    out: Dict[str, int] = {}
    for c in gtf_chroms:
        for cand in (c, "chr" + c, c.replace("chr", ""),
                     "chrM" if c in ("MT", "M") else None,
                     "MT" if c in ("chrM", "M") else None):
            if cand is not None and cand in idx:
                out[c] = idx[cand]
                break
    return out


class Annotation:
    """Conteneur des tableaux d'annotation en coordonnees globales."""

    FIELDS = ("site_pos", "site_kind", "site_gene", "intron_start", "intron_end",
              "intron_gene", "exon_start", "exon_end", "exon_gene",
              "gene_start", "gene_end", "gene_strand", "gene_chrom")

    def __init__(self, **kw):
        for f in self.FIELDS:
            setattr(self, f, kw[f])
        self.gene_name: List[str] = list(kw["gene_name"])
        self.gene_id: List[str] = list(kw["gene_id"])

    @property
    def n_sites(self) -> int:
        return len(self.site_pos)

    def save(self, path: Path) -> None:
        d = {f: getattr(self, f) for f in self.FIELDS}
        d["gene_name"] = np.array(self.gene_name, dtype=object)
        d["gene_id"] = np.array(self.gene_id, dtype=object)
        np.savez_compressed(winpath(path), **d)

    @classmethod
    def load(cls, path: Path) -> "Annotation":
        z = np.load(winpath(path), allow_pickle=True)
        kw = {f: z[f] for f in cls.FIELDS}
        kw["gene_name"] = z["gene_name"].tolist()
        kw["gene_id"] = z["gene_id"].tolist()
        return cls(**kw)

    @classmethod
    def empty(cls) -> "Annotation":
        z64 = np.zeros(0, dtype=np.int64)
        return cls(site_pos=z64, site_kind=np.zeros(0, np.int8),
                   site_gene=np.zeros(0, np.int32),
                   intron_start=z64, intron_end=z64.copy(),
                   intron_gene=np.zeros(0, np.int32),
                   exon_start=z64.copy(), exon_end=z64.copy(),
                   exon_gene=np.zeros(0, np.int32),
                   gene_start=z64.copy(), gene_end=z64.copy(),
                   gene_strand=np.zeros(0, np.int8),
                   gene_chrom=np.zeros(0, np.int32),
                   gene_name=[], gene_id=[])


def parse_gtf(gtf: str, ref_names: Sequence[str], ref_off: np.ndarray,
              feature_exon: str = "exon",
              keep_biotypes: Optional[Sequence[str]] = None) -> Annotation:
    log(f"Lecture du GTF : {gtf}")
    gene_key: Dict[str, int] = {}
    gene_name: List[str] = []
    gene_id: List[str] = []
    g_chrom: List[int] = []
    g_start: List[int] = []
    g_end: List[int] = []
    g_strand: List[int] = []

    tx_key: Dict[str, int] = {}
    ex_tx = array("i")
    ex_s = array("q")
    ex_e = array("q")
    ex_g = array("i")

    unknown_chrom: Dict[str, int] = {}
    idx_by_name = {n: i for i, n in enumerate(ref_names)}
    chrom_cache: Dict[str, int] = {}
    n_lines = 0

    with open_maybe_gzip(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            n_lines += 1
            f = line.split("\t")
            if len(f) < 9:
                continue
            feat = f[2]
            if feat != feature_exon and feat != "gene":
                continue
            c = f[0]
            ri = chrom_cache.get(c, -2)
            if ri == -2:
                ri = -1
                for cand in (c, "chr" + c, c[3:] if c.startswith("chr") else None,
                             "chrM" if c in ("MT", "M") else None,
                             "MT" if c == "chrM" else None):
                    if cand and cand in idx_by_name:
                        ri = idx_by_name[cand]
                        break
                chrom_cache[c] = ri
            if ri < 0:
                unknown_chrom[c] = unknown_chrom.get(c, 0) + 1
                continue
            attrs = f[8]
            gid = _attr(attrs, "gene_id") or "NA"
            gi = gene_key.get(gid, -1)
            if gi < 0:
                gi = len(gene_id)
                gene_key[gid] = gi
                gene_id.append(gid)
                gene_name.append(_attr(attrs, "gene_name") or gid)
                g_chrom.append(ri)
                g_start.append(1 << 62)
                g_end.append(-1)
                g_strand.append(1 if f[6] == "+" else (-1 if f[6] == "-" else 0))
            gs = int(f[3]) - 1 + int(ref_off[ri])   # 0-based global
            ge = int(f[4]) + int(ref_off[ri])       # exclusif
            if gs < g_start[gi]:
                g_start[gi] = gs
            if ge > g_end[gi]:
                g_end[gi] = ge
            if feat != feature_exon:
                continue
            if keep_biotypes:
                bt = _attr(attrs, "gene_type") or _attr(attrs, "gene_biotype")
                if bt is not None and bt not in keep_biotypes:
                    continue
            tid = _attr(attrs, "transcript_id") or gid
            ti = tx_key.get(tid, -1)
            if ti < 0:
                ti = len(tx_key)
                tx_key[tid] = ti
            ex_tx.append(ti)
            ex_s.append(gs)
            ex_e.append(ge)
            ex_g.append(gi)

    if unknown_chrom:
        top = sorted(unknown_chrom.items(), key=lambda kv: -kv[1])[:5]
        warn(f"{len(unknown_chrom)} chromosomes du GTF absents du BAM (ignores), "
             f"ex: {[t[0] for t in top]}")
    log(f"GTF : {n_lines} lignes utiles, {len(gene_id)} genes, "
        f"{len(ex_s)} exons, {len(tx_key)} transcrits")

    ex_tx_a = np.frombuffer(ex_tx, dtype=np.int32).copy()
    ex_s_a = np.frombuffer(ex_s, dtype=np.int64).copy()
    ex_e_a = np.frombuffer(ex_e, dtype=np.int64).copy()
    ex_g_a = np.frombuffer(ex_g, dtype=np.int32).copy()
    del ex_tx, ex_s, ex_e, ex_g

    if len(ex_s_a) == 0:
        warn("Aucun exon lu dans le GTF : verifier le format / les chromosomes.")
        return Annotation.empty()

    # introns = intervalles entre exons consecutifs d'un meme transcrit
    order = np.lexsort((ex_s_a, ex_tx_a))
    t = ex_tx_a[order]
    s = ex_s_a[order]
    e = ex_e_a[order]
    g = ex_g_a[order]
    same = t[1:] == t[:-1]
    istart = e[:-1]
    iend = s[1:]
    ok = same & (iend > istart)
    intron_start = istart[ok]
    intron_end = iend[ok]
    intron_gene = g[:-1][ok]

    # dedupliquer les introns
    key = intron_start * np.int64(1 << 21) + np.minimum(
        intron_end - intron_start, (1 << 21) - 1)
    uk, ui = np.unique(key, return_index=True)
    intron_start = intron_start[ui]
    intron_end = intron_end[ui]
    intron_gene = intron_gene[ui]

    # sites d'epissage = bornes des introns
    allpos = np.concatenate([intron_start, intron_end])
    allkind = np.concatenate([np.full(len(intron_start), 1, np.int8),
                              np.full(len(intron_end), 2, np.int8)])
    allgene = np.concatenate([intron_gene, intron_gene])
    o = np.argsort(allpos, kind="stable")
    allpos = allpos[o]
    allkind = allkind[o]
    allgene = allgene[o]
    uniq, first = np.unique(allpos, return_index=True)
    # fusionner le type (donneur/accepteur) pour les positions dupliquees
    kind = np.zeros(len(uniq), dtype=np.int8)
    np.bitwise_or.at(kind, np.searchsorted(uniq, allpos), allkind)
    site_pos = uniq
    site_kind = kind
    site_gene = allgene[first].astype(np.int32)

    # exons uniques (pour le saut d'exon)
    ekey = ex_s_a * np.int64(1 << 21) + np.minimum(ex_e_a - ex_s_a, (1 << 21) - 1)
    _, eu = np.unique(ekey, return_index=True)
    exon_start = ex_s_a[eu]
    exon_end = ex_e_a[eu]
    exon_gene = ex_g_a[eu]
    eo = np.argsort(exon_start, kind="stable")

    log(f"Annotation : {len(site_pos)} sites d'epissage, "
        f"{len(intron_start)} introns, {len(exon_start)} exons uniques")

    return Annotation(
        site_pos=site_pos.astype(np.int64),
        site_kind=site_kind,
        site_gene=site_gene,
        intron_start=intron_start.astype(np.int64),
        intron_end=intron_end.astype(np.int64),
        intron_gene=intron_gene.astype(np.int32),
        exon_start=exon_start[eo].astype(np.int64),
        exon_end=exon_end[eo].astype(np.int64),
        exon_gene=exon_gene[eo].astype(np.int32),
        gene_start=np.array(g_start, dtype=np.int64),
        gene_end=np.array(g_end, dtype=np.int64),
        gene_strand=np.array(g_strand, dtype=np.int8),
        gene_chrom=np.array(g_chrom, dtype=np.int32),
        gene_name=gene_name,
        gene_id=gene_id,
    )


def canonical_introns(gtf: str, ref_names: Sequence[str], ref_off: np.ndarray,
                      cache: Optional[Path] = None, tag: str = "Ensembl_canonical"):
    """Introns des transcrits marques `Ensembl_canonical` (GENCODE / Ensembl).

    Renvoie un tableau trie de cles (debut global << 24 | longueur), meme
    convention que les jonctions. Mis en cache (le GTF n'est relu que s'il
    change). Si le GTF ne contient aucun tag `Ensembl_canonical`, on se rabat
    sur `MANE_Select`, puis sur aucun (tableau vide).
    """
    import json as _json
    import os as _os
    sig = ""
    try:
        st = _os.stat(gtf)
        sig = f"{_os.path.abspath(gtf)}|{st.st_size}|{int(st.st_mtime)}|{tag}"
    except OSError:
        return None
    if cache is not None and Path(cache).exists():
        try:
            z = np.load(winpath(cache), allow_pickle=False)
            if str(z["sig"]) == sig:
                return z["keys"]
        except Exception:  # noqa: BLE001
            pass
    log(f"Introns des transcrits {tag} : lecture de {gtf}")
    idx_by_name = {n: i for i, n in enumerate(ref_names)}
    tx: Dict[str, list] = {}
    mane: Dict[str, list] = {}
    with open_maybe_gzip(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon":
                continue
            a = f[8]
            is_c = f'tag "{tag}"' in a
            is_m = 'tag "MANE_Select"' in a
            if not (is_c or is_m):
                continue
            c = f[0]
            ri = idx_by_name.get(c, idx_by_name.get("chr" + c, idx_by_name.get(
                c[3:] if c.startswith("chr") else "", -1)))
            if ri < 0:
                continue
            tid = _attr(a, "transcript_id") or ""
            ex = (int(f[3]) - 1 + int(ref_off[ri]), int(f[4]) + int(ref_off[ri]))
            if is_c:
                tx.setdefault(tid, []).append(ex)
            if is_m:
                mane.setdefault(tid, []).append(ex)
    src = tx if tx else mane
    if not tx and mane:
        warn(f"aucun tag {tag} dans le GTF : transcrits MANE_Select utilises")
    keys = []
    for exs in src.values():
        exs.sort()
        for (s1, e1), (s2, _e2) in zip(exs[:-1], exs[1:]):
            if s2 > e1:
                keys.append((e1 << 24) | ((s2 - e1) & ((1 << 24) - 1)))
    arr = np.unique(np.array(keys, dtype=np.int64))
    log(f"{len(src)} transcrits canoniques, {len(arr)} introns")
    if cache is not None:
        try:
            with open(winpath(cache), "wb") as fh:
                np.savez(fh, keys=arr, sig=np.array(sig))
        except OSError:
            pass
    return arr
