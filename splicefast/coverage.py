"""
Couverture des introns : mesure rapide, sur TOUS les introns annotes.

Principe (equivalent mathematique de PanDepth, mais faisable en Python pur)
--------------------------------------------------------------------------
PanDepth alloue un tableau de profondeur par chromosome et incremente chaque
base couverte par chaque read. En C c'est imbattable ; en Python ce serait des
milliards d'increments. On calcule donc exactement les memes quantites sans
jamais materialiser la profondeur base par base :

  * les points d'interet (bornes d'introns et d'exons) definissent des
    SEGMENTS. Pour chaque bloc aligne [a, b) on ajoute en O(1) :
        - la longueur partielle dans le segment de a et dans celui de b,
        - +1 / -1 dans un tableau de differences pour les segments entierement
          couverts (leur contribution = largeur du segment).
    Une seule somme cumulee a la fin donne le TOTAL DE BASES par segment,
    donc, par sommes prefixes, la somme de profondeur sur n'importe quel
    intervalle.
  * la profondeur en un POINT precis (la derniere base de l'exon amont, la
    premiere de l'exon aval) vient d'un second tableau de differences :
    +1 en a, -1 en b, somme cumulee.

Tout est vectorise par lots de ~500 000 blocs : le cout par read redevient
celui du decodage BAM, pas celui de numpy appele un million de fois.

Normalisation
-------------
    ratio = profondeur moyenne de l'intron / profondeur moyenne des DEUX
            bases qui le bordent (derniere base de l'exon amont, premiere
            base de l'exon aval)

On ne prend PAS la moyenne de l'exon entier : une bosse de couverture au
milieu d'un exon (3'UTR, exon alternatif, artefact d'amorcage) ecraserait le
signal. Les deux bases bordantes sont la reference la plus locale possible.

Les positions exoniques d'AUTRES transcrits (exons alternatifs, genes
imbriques) sont exclues du corps de l'intron : sans cela, un exon interne se
lirait comme une retention.
"""
from __future__ import annotations

import hashlib
import os
import struct
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .bam import DEFAULT_EXCLUDE, ref_offsets
from .bgzf import BgzfBlockStream, split_bam, vsplit
from .util import Progress, log, warn, winpath

_CORE = struct.Struct("<iiBBHHHi")
_U32 = struct.Struct("<I")
FLUSH = 400_000          # blocs alignes accumules avant vectorisation


# --------------------------------------------------------------------------- #
# Cibles
# --------------------------------------------------------------------------- #
def _merge_intervals(s: np.ndarray, e: np.ndarray):
    """Fusionne des intervalles [s, e) qui se chevauchent ou se touchent
    (vectorise : plus de boucle Python sur les 600 000 exons de GENCODE)."""
    if len(s) == 0:
        return s, e
    o = np.argsort(s, kind="stable")
    s, e = s[o], e[o]
    run_e = np.maximum.accumulate(e)
    new = np.ones(len(s), dtype=bool)
    new[1:] = s[1:] > run_e[:-1]
    gid = np.cumsum(new) - 1
    ms = s[new]
    me = np.zeros(len(ms), dtype=np.int64)
    np.maximum.at(me, gid, e)
    return ms.astype(np.int64), me.astype(np.int64)


class CovTargets:
    """Introns a mesurer + grille de segments partagee."""

    def __init__(self, q, seg_w, seg_free, bpos, chrom, start, end, gene,
                 qi_s, qi_e, bi_l, bi_r):
        self.q = q                  # bornes des segments (global, trie, unique)
        self.seg_w = seg_w          # largeur de chaque segment
        self.seg_free = seg_free    # largeur NON exonique de chaque segment
        self.bpos = bpos            # positions bordantes (trie, unique)
        self.chrom = chrom
        self.start = start          # coordonnees locales 0-based
        self.end = end
        self.gene = gene
        self.qi_s = qi_s            # index du segment ou commence l'intron
        self.qi_e = qi_e
        self.bi_l = bi_l            # index dans bpos de la base bordante gauche
        self.bi_r = bi_r

    @property
    def n(self) -> int:
        return len(self.start)


def build_cov_targets(ann, ref_names: Sequence[str], ref_lengths: Sequence[int],
                      genes: Optional[Sequence[str]] = None,
                      edge: int = 10, min_len: int = 60,
                      max_len: int = 1_000_000,
                      border_width: int = 25) -> CovTargets:
    ref_off = ref_offsets(ref_lengths)
    istart = np.asarray(ann.intron_start, dtype=np.int64)
    iend = np.asarray(ann.intron_end, dtype=np.int64)
    igene = np.asarray(ann.intron_gene, dtype=np.int64)

    keep = (iend - istart >= min_len) & (iend - istart <= max_len)
    if genes:
        want = {g.strip().upper() for g in genes if g.strip()}
        gsel = np.array([str(g).upper() in want for g in ann.gene_name])
        keep &= gsel[np.clip(igene, 0, len(gsel) - 1)] & (igene >= 0)
        log(f"couverture : restriction a {int(gsel.sum())} gene(s)")
    # un intron peut etre decrit par plusieurs transcrits : une seule mesure
    idx = np.where(keep)[0]
    uniq, first = np.unique(istart[idx] * (1 << 32) + (iend[idx] & 0xFFFFFFFF),
                            return_index=True)
    idx = idx[np.sort(first)]
    istart, iend, igene = istart[idx], iend[idx], igene[idx]

    ex_s, ex_e = _merge_intervals(np.asarray(ann.exon_start, dtype=np.int64),
                                  np.asarray(ann.exon_end, dtype=np.int64))

    # corps mesure de l'intron : [start+edge, end-edge)
    cs = istart + edge
    ce = iend - edge
    good = ce - cs >= 20
    istart, iend, igene, cs, ce = (istart[good], iend[good], igene[good],
                                   cs[good], ce[good])

    # grille de segments : bornes des corps d'introns ET des exons
    W = max(1, int(border_width))
    # istart-1 et iend+1 isolent les bases bordantes dans des segments de
    # largeur 1 : leur profondeur = bases du segment (utilise par sfscan).
    # Raffiner la grille ne change aucune somme.
    q = np.unique(np.concatenate([cs, ce, ex_s, ex_e,
                                  istart - W, istart - 1, istart, iend,
                                  iend + 1, iend + W]))
    seg_w = np.diff(q)
    # largeur non exonique de chaque segment
    mid = q[:-1]
    i = np.searchsorted(ex_s, mid, side="right") - 1
    inside = (i >= 0) & (mid < ex_e[np.clip(i, 0, len(ex_e) - 1)])
    seg_free = np.where(inside, 0, seg_w).astype(np.int64)

    bpos = np.unique(np.concatenate([istart - 1, iend]))

    tg = CovTargets(
        q=q, seg_w=seg_w.astype(np.int64), seg_free=seg_free, bpos=bpos,
        chrom=None, start=None, end=None, gene=igene,
        qi_s=np.searchsorted(q, cs, side="left"),
        qi_e=np.searchsorted(q, ce, side="left"),
        bi_l=np.searchsorted(bpos, istart - 1),
        bi_r=np.searchsorted(bpos, iend))
    # fenetre exonique bordante : W bases de part et d'autre de l'intron
    tg.wi_l0 = np.searchsorted(q, istart - W, side="left")
    tg.wi_l1 = np.searchsorted(q, istart, side="left")
    tg.wi_r0 = np.searchsorted(q, iend, side="left")
    tg.wi_r1 = np.searchsorted(q, iend + W, side="left")
    tg.border_width = W
    tg.pt_l = np.searchsorted(q, istart - 1, side="left")
    tg.pt_r = np.searchsorted(q, iend, side="left")
    tg.idx = np.ascontiguousarray(np.stack(
        [tg.qi_s, tg.qi_e, tg.wi_l0, tg.wi_l1, tg.wi_r0, tg.wi_r1,
         tg.pt_l, tg.pt_r], axis=1).astype(np.int64))
    h = hashlib.md5()
    for arr in (q, seg_free, tg.idx):
        h.update(np.ascontiguousarray(arr, dtype=np.int64).tobytes())
    tg.hash = h.hexdigest()[:16]
    # coordonnees locales pour l'affichage
    gc = np.asarray(ann.gene_chrom, dtype=np.int64)
    refid = gc[np.clip(igene, 0, len(gc) - 1)]
    off = np.asarray(ref_off, dtype=np.int64)[np.clip(refid, 0, len(ref_off) - 1)]
    tg.chrom = np.array([ref_names[r] if 0 <= r < len(ref_names) else "?"
                         for r in refid], dtype=object)
    tg.start = istart - off
    tg.end = iend - off
    log(f"couverture : {tg.n} introns, {len(seg_w)} segments, "
        f"{len(bpos)} bases bordantes")
    return tg


# --------------------------------------------------------------------------- #
# Balayage d'une tranche de BAM
# --------------------------------------------------------------------------- #
def scan_chunk_cov(bam: str, vstart: int, vend: int, ref_off: np.ndarray,
                   q: np.ndarray, bpos: np.ndarray, min_mapq: int = 0,
                   exclude_flags: int = DEFAULT_EXCLUDE,
                   read_buffer: int = 4 << 20):
    """Bases par segment et profondeur aux bases bordantes, pour [vstart, vend)."""
    nseg = len(q) - 1
    part = np.zeros(nseg + 2, dtype=np.int64)     # contributions partielles
    full = np.zeros(nseg + 2, dtype=np.int64)     # segments entierement couverts
    dep = np.zeros(len(bpos) + 2, dtype=np.int64)
    A: List[int] = []
    B: List[int] = []
    n_reads = 0

    def flush():
        if not A:
            return
        a = np.fromiter(A, dtype=np.int64, count=len(A))
        b = np.fromiter(B, dtype=np.int64, count=len(B))
        A.clear()
        B.clear()
        # --- profondeur aux bases bordantes : a <= x < b
        i0 = np.searchsorted(bpos, a, side="left")
        i1 = np.searchsorted(bpos, b, side="left")
        m = i1 > i0
        if m.any():
            np.add.at(dep, i0[m], 1)
            np.add.at(dep, i1[m], -1)
        # --- bases par segment
        ia = np.searchsorted(q, a, side="right") - 1
        ib = np.searchsorted(q, b, side="right") - 1
        same = (ia == ib) & (ia >= 0) & (ia < nseg)
        if same.any():
            np.add.at(part, ia[same], (b - a)[same])
        diff = ~(ia == ib)
        if diff.any():
            va = diff & (ia >= 0) & (ia < nseg)
            if va.any():
                np.add.at(part, ia[va], (q[np.clip(ia[va] + 1, 0, nseg)]
                                         - a[va]))
            vb = diff & (ib >= 0) & (ib < nseg)
            if vb.any():
                np.add.at(part, ib[vb], (b[vb] - q[np.clip(ib[vb], 0, nseg)]))
            lo = np.clip(ia[diff] + 1, 0, nseg)
            hi = np.clip(ib[diff], 0, nseg)
            k = hi > lo
            if k.any():
                np.add.at(full, lo[k], 1)
                np.add.at(full, hi[k], -1)

    start_c, start_u = vsplit(vstart)
    stop_c, stop_u = vsplit(vend) if vend else (None, None)
    buf = bytearray()
    stream = BgzfBlockStream(bam, start_c, read_buffer=read_buffer)
    first = True
    finished = False
    try:
        for coff, data in stream.blocks():
            if first:
                buf += data[start_u:]
                first = False
            else:
                buf += data
            if stop_c is not None and coff >= stop_c:
                finished = True
            off = 0
            while True:
                if len(buf) - off < 4:
                    break
                (bs,) = struct.unpack_from("<i", buf, off)
                if len(buf) - off < bs + 4:
                    break
                (refid, pos, l_rn, mapq, _bin, n_cig, flag,
                 _l_seq) = _CORE.unpack_from(buf, off + 4)
                off_next = off + 4 + bs
                if refid < 0 or (flag & exclude_flags) or mapq < min_mapq \
                        or n_cig == 0:
                    off = off_next
                    continue
                n_reads += 1
                cur = int(ref_off[refid]) + pos
                cig_off = off + 36 + l_rn
                for i in range(n_cig):
                    (c,) = _U32.unpack_from(buf, cig_off + 4 * i)
                    op = c & 0xF
                    ln = c >> 4
                    if op == 0 or op == 7 or op == 8:
                        A.append(cur)
                        cur += ln
                        B.append(cur)
                    elif op == 2 or op == 3:
                        cur += ln
                if len(A) >= FLUSH:
                    flush()
                off = off_next
            if off:
                del buf[:off]
            if finished:
                break
    finally:
        stream.close()
    flush()
    base = part[:nseg] + np.cumsum(full[:nseg]) * (q[1:] - q[:-1])
    return base, np.cumsum(dep[:len(bpos)]), n_reads


def _task(args):
    (bam, vstart, vend, ref_off, qpath, bpath, mapq, excl) = args
    q = np.load(qpath, mmap_mode="r")
    bpos = np.load(bpath, mmap_mode="r")
    return scan_chunk_cov(bam, vstart, vend, np.asarray(ref_off),
                          np.asarray(q), np.asarray(bpos),
                          min_mapq=mapq, exclude_flags=excl)


def _finish_python(tg: CovTargets, base: np.ndarray, dep: np.ndarray) -> np.ndarray:
    """Profondeurs [intron, base bordante, fenetre bordante] (M x 3)."""
    cum_w = np.concatenate([[0], np.cumsum(tg.seg_free)])
    ex_w = tg.seg_w - tg.seg_free                  # largeur exonique
    cum_ex_w = np.concatenate([[0], np.cumsum(ex_w)])
    # somme de profondeur sur la partie NON exonique du corps de l'intron
    free_base = np.where(tg.seg_free > 0, base, 0)
    cum_b = np.concatenate([[0], np.cumsum(free_base)])
    s_sum = cum_b[tg.qi_e] - cum_b[tg.qi_s]
    s_len = cum_w[tg.qi_e] - cum_w[tg.qi_s]
    out = np.empty((tg.n, 3), dtype=np.float32)
    out[:, 0] = s_sum / np.maximum(s_len, 1)
    out[:, 1] = (dep[tg.bi_l] + dep[tg.bi_r]) / 2.0
    # profondeur moyenne sur les W bases EXONIQUES qui bordent l'intron
    cum_eb = np.concatenate([[0], np.cumsum(np.where(tg.seg_free > 0, 0, base))])
    wb = ((cum_eb[tg.wi_l1] - cum_eb[tg.wi_l0])
          + (cum_eb[tg.wi_r1] - cum_eb[tg.wi_r0]))
    ww = ((cum_ex_w[tg.wi_l1] - cum_ex_w[tg.wi_l0])
          + (cum_ex_w[tg.wi_r1] - cum_ex_w[tg.wi_r0]))
    out[:, 2] = wb / np.maximum(ww, 1)
    return out


def sample_coverage_python(wd, bam: str, tg: CovTargets, ref_off: np.ndarray,
                           threads: int = 4, min_mapq: int = 0,
                           exclude_flags: int = DEFAULT_EXCLUDE,
                           chunk_mb: int = 64) -> np.ndarray:
    """Ancien moteur, en Python pur (repli si sfscan n'est pas compile)."""
    qpath = winpath(wd.tmp / "cov_q.npy")
    bpath = winpath(wd.tmp / "cov_b.npy")
    wd.tmp.mkdir(parents=True, exist_ok=True)
    np.save(qpath, tg.q)
    np.save(bpath, tg.bpos)
    chunks = split_bam(bam, max(1, threads), chunk_mb << 20)
    args = [(bam, v0, v1, ref_off, qpath, bpath, min_mapq, exclude_flags)
            for v0, v1 in chunks]
    if threads <= 1 or len(args) == 1:
        res = [_task(a) for a in args]
    else:
        res = []
        with ProcessPoolExecutor(max_workers=threads) as ex:
            futs = [ex.submit(_task, a) for a in args]
            for f in as_completed(futs):
                res.append(f.result())
    base = dep = None
    for b, d, _n in res:
        base = b if base is None else base + b
        dep = d if dep is None else dep + d
    return _finish_python(tg, base, dep)


# --------------------------------------------------------------------------- #
# Cache par echantillon : un run interrompu reprend la ou il s'est arrete
# --------------------------------------------------------------------------- #
def cov_sample_hash(tg: CovTargets, bam: str, min_mapq: int,
                    exclude_flags: int) -> str:
    h = hashlib.md5(f"{tg.hash}|{bam}|{min_mapq}|{exclude_flags}".encode())
    return h.hexdigest()[:16]


def cov_cache_path(wd, sid: str):
    from .util import safe_name
    return wd.covcounts / f"{safe_name(sid)}.npz"


def cov_cache_ok(wd, sid: str, hsh: str, M: int) -> bool:
    p = cov_cache_path(wd, sid)
    if not p.exists():
        return False
    try:
        with np.load(winpath(p), allow_pickle=False) as z:
            return str(z["hash"]) == hsh and z["res"].shape == (M, 3)
    except Exception:
        return False


def cov_cache_save(wd, sid: str, hsh: str, res: np.ndarray) -> None:
    p = cov_cache_path(wd, sid)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".part")
    with open(winpath(tmp), "wb") as fh:
        np.savez(fh, res=np.asarray(res, dtype=np.float32), hash=np.array(hsh))
    os.replace(winpath(tmp), winpath(p))


def write_cov_grid(wd, tg: CovTargets) -> Dict[str, str]:
    """Grille de segments pour sfscan (ecrite une fois par jeu de cibles)."""
    d = wd.tmp / f"covgrid_{tg.hash}"
    d.mkdir(parents=True, exist_ok=True)
    paths = {"grid": d / "q.npy", "free": d / "free.npy", "idx": d / "idx.npy"}
    if not all(p.exists() for p in paths.values()):
        np.save(paths["grid"], np.asarray(tg.q, dtype=np.int64))
        np.save(paths["free"], np.asarray(tg.seg_free, dtype=np.int64))
        np.save(paths["idx"], tg.idx)
    return {k: str(v) for k, v in paths.items()}


def cov_scan_args(wd, tg: CovTargets, sid: str, min_mapq: int,
                  exclude_flags: int) -> Tuple[List[str], str]:
    g = write_cov_grid(wd, tg)
    from .util import safe_name
    out = wd.tmp / f"cov_{safe_name(sid)}.npy"
    return (["--cov", str(out), "--cov-grid", g["grid"], "--cov-free", g["free"],
             "--cov-idx", g["idx"], "--cov-mapq", str(min_mapq),
             "--cov-exclude", str(exclude_flags)], str(out))


def measure_coverage(wd, samples, tg: CovTargets, ref_lengths,
                     threads: int = 4, min_mapq: int = 0,
                     exclude_flags: int = DEFAULT_EXCLUDE,
                     chunk_mb: int = 64, jobs: int = 1, engine: str = "auto",
                     ref_names: Optional[Sequence[str]] = None):
    """Profondeur moyenne intronique et bordante, pour tous les echantillons.

    Chaque echantillon mesure est mis en cache (workdir/covcounts) : relancer
    la commande apres une coupure ne relit que les BAM manquants. Les
    echantillons deja mesures pendant `count --retention` sont reutilises.
    """
    from concurrent.futures import ThreadPoolExecutor
    from .scanner import (ScanError, resolve_engine, run_sfscan, split_threads,
                          write_ref_files)
    from .util import human

    ref_off = ref_offsets(ref_lengths)
    S = len(samples)
    M = tg.n
    hashes = {sid: cov_sample_hash(tg, bam, min_mapq, exclude_flags)
              for sid, bam in samples}
    todo = [(sid, bam) for sid, bam in samples
            if not cov_cache_ok(wd, sid, hashes[sid], M)]
    log(f"couverture : {len(todo)} echantillon(s) a lire, "
        f"{S - len(todo)} deja mesure(s)")
    failed: List[str] = []
    if todo:
        eng = resolve_engine(engine)
        tot = sum(os.path.getsize(b) for _, b in todo if os.path.exists(b))
        log(f"couverture : moteur {eng}, {human(tot)} de BAM")
        pr = Progress(len(todo), "couverture des introns")
        if eng == "cpp":
            jobs, per = split_threads(threads, jobs)
            refs = write_ref_files(wd.tmp, ref_names or [], ref_off)

            def one(item):
                sid, bam = item
                extra, out = cov_scan_args(wd, tg, sid, min_mapq, exclude_flags)
                args = (["--bam", bam, "--ref-off", refs["ref_off"], "-@", per]
                        + (["--ref-names", refs["ref_names"]] if ref_names else [])
                        + extra)
                info = run_sfscan(args, wait_path=bam, label=sid)
                res = np.load(out)
                os.remove(out)
                cov_cache_save(wd, sid, hashes[sid], res)
                return sid, info

            with ThreadPoolExecutor(max_workers=jobs) as ex:
                futs = {ex.submit(one, it): it[0] for it in todo}
                from concurrent.futures import as_completed as _ac
                for f in _ac(futs):
                    sid = futs[f]
                    try:
                        _, info = f.result()
                        log(f"  {sid} : {info.get('n_cov_reads', 0):,} reads en "
                            f"{info.get('seconds', 0):.0f} s")
                    except ScanError as e:
                        warn(f"{sid} : ECHEC couverture ({e}) - mis de cote")
                        failed.append(sid)
                    pr.update()
        else:
            from .scanner import header_matches
            for sid, bam in todo:
                if not header_matches(bam, ref_names):
                    warn(f"{sid} : chromosomes differents du projet, "
                         f"illisible par le lecteur Python (compiler sfscan)")
                    failed.append(sid)
                    pr.update()
                    continue
                try:
                    res = sample_coverage_python(wd, bam, tg, ref_off, threads,
                                                 min_mapq, exclude_flags,
                                                 chunk_mb)
                    cov_cache_save(wd, sid, hashes[sid], res)
                except Exception as e:  # noqa: BLE001
                    warn(f"{sid} : ECHEC couverture ({e.__class__.__name__}: "
                         f"{e}) - mis de cote")
                    failed.append(sid)
                pr.update()
        pr.close()
    if failed:
        raise SystemExit(f"{len(failed)} echantillon(s) sans couverture : "
                         f"{', '.join(failed[:8])}. Relancer la meme commande "
                         f"pour les reprendre (les autres sont gardes).")

    # assemblage : memmap plutot que RAM (147 x 1,2 M introns = 700 Mo/matrice)
    D = np.memmap(winpath(wd.store / "cov.D.dat"), dtype=np.float32,
                  mode="w+", shape=(S, M))
    Bd = np.memmap(winpath(wd.store / "cov.Bd.dat"), dtype=np.float32,
                   mode="w+", shape=(S, M))
    Bw = np.memmap(winpath(wd.store / "cov.Bw.dat"), dtype=np.float32,
                   mode="w+", shape=(S, M))
    for si, (sid, _bam) in enumerate(samples):
        with np.load(winpath(cov_cache_path(wd, sid))) as z:
            r = z["res"]
            D[si], Bd[si], Bw[si] = r[:, 0], r[:, 1], r[:, 2]
    for arr in (D, Bd, Bw):
        arr.flush()
    del D, Bd, Bw
    return S, M


def save_coverage(wd, tg: CovTargets, samples, S: int, M: int) -> None:
    path = wd.store / "coverage.npz"
    with open(winpath(path), "wb") as fh:
        np.savez_compressed(
            fh, chrom=np.array(list(tg.chrom), dtype=object),
            start=tg.start, end=tg.end, gene=tg.gene,
            length=(tg.end - tg.start), shape=np.array([S, M]),
            border_width=np.array(getattr(tg, "border_width", 1)),
            samples=np.array([s for s, _ in samples], dtype=object))
    log(f"couverture enregistree : {path} "
        f"({S} echantillons x {M} introns)")


def load_coverage(wd):
    """Coordonnees + acces memmap aux matrices de profondeur."""
    path = wd.store / "coverage.npz"
    if not path.exists():
        return None
    z = np.load(winpath(path), allow_pickle=True)
    d = {k: z[k] for k in z.files}
    S, M = (int(x) for x in d["shape"])
    d["D"] = np.memmap(winpath(wd.store / "cov.D.dat"), dtype=np.float32,
                       mode="r", shape=(S, M))
    d["Bd"] = np.memmap(winpath(wd.store / "cov.Bd.dat"), dtype=np.float32,
                        mode="r", shape=(S, M))
    pw = wd.store / "cov.Bw.dat"
    d["Bw"] = (np.memmap(winpath(pw), dtype=np.float32, mode="r",
                         shape=(S, M)) if pw.exists() else d["Bd"])
    return d


# --------------------------------------------------------------------------- #
# Test statistique et sortie
# --------------------------------------------------------------------------- #
def coverage_model(wd, min_border: float = 10.0, rho_floor: float = 0.005,
                   min_samples: int = 25,
                   n_factors: Optional[int] = None) -> Optional[dict]:
    """Teste le rapport profondeur(intron)/profondeur(bases bordantes).

    Le rapport est converti en proportion  psi = D / (D + Bd)  pour entrer
    dans le test beta-binomial : psi = 0.5 signifie un intron aussi couvert
    que ses bords (retention complete), psi = 0.33 environ 50 % de retention.
    K et N sont des profondeurs, donc des NOMBRES DE READS (la profondeur en
    une base = nombre de reads qui la couvrent) : pas de sur-estimation de
    l'information comme lorsqu'on comptait les bases.
    """
    from .stats import (bb_logp_screened, bh, design_matrix, latent_factors,
                        refine_rho, rho_moments, robust_expected, trim_mask)

    from .store import read_reference

    d = load_coverage(wd)
    if d is None:
        return None
    Dm, Bm = d["D"], d["Bw"]
    samples = [str(s) for s in d["samples"]]
    S, M = Dm.shape
    ref = read_reference(wd, samples)          # qui definit la norme
    if ref.sum() < 5:
        ref = np.ones(S, dtype=bool)
    if not ref.all():
        log(f"couverture : norme definie par {int(ref.sum())} echantillon(s), "
            f"{int((~ref).sum())} teste(s) hors reference")

    # selection des introns mesurables, par blocs (la matrice reste sur disque)
    # Un intron est testable des lors qu'assez d'echantillons y ont une
    # couverture exploitable -- PAS quand la mediane de toute la cohorte est
    # elevee. Un gene a expression tres variable (SMC1A : parfois 5 reads,
    # parfois 200) serait sinon ecarte alors que 40 echantillons suffisent
    # largement a definir la norme.
    # sur une petite cohorte, on abaisse l'exigence plutot que de tout exclure
    n_ref = int(ref.sum())
    need = int(min(min_samples, max(5, n_ref // 3)))
    cstep = max(1000, 8_000_000 // max(S, 1))
    keep_cols = []
    for a in range(0, M, cstep):
        b = min(a + cstep, M)
        n_ok = (np.asarray(Bm[:, a:b])[ref] >= min_border).sum(axis=0)
        keep_cols.append(np.where(n_ok >= need)[0] + a)
    cols = np.concatenate(keep_cols) if keep_cols else np.zeros(0, int)
    if not len(cols):
        warn(f"couverture : aucun intron avec >= {need} echantillons "
             f"couverts (bordure >= {min_border})")
        return None
    M2 = len(cols)

    q = int(n_factors) if n_factors is not None else (
        max(1, min(4, n_ref // 12)) if M2 >= 200 else 0)
    q = min(q, max(0, n_ref - 3))
    if q > 0:
        sub = np.arange(M2)
        if M2 > 20000:
            sub = np.sort(np.random.default_rng(0).choice(M2, 20000,
                                                          replace=False))
        cs = np.sort(cols[sub])
        kk = np.asarray(Dm[:, cs], dtype=np.float64)
        bb = np.asarray(Bm[:, cs], dtype=np.float64)
        # seuls les introns couverts PARTOUT servent aux facteurs latents
        full = (bb[ref] >= min_border).all(axis=0)
        if full.sum() >= 50:
            kk, bb = kk[:, full], bb[:, full]
        Z, _e = latent_factors(kk / np.maximum(kk + bb, 1e-9), q, ref=ref)
        del kk, bb
        X = design_matrix(Z[:, :q])
        log(f"couverture : {q} facteur(s) latent(s), {M2} introns testables")
    else:
        X = design_matrix(np.zeros((S, 0)))
        log(f"couverture : modele sans facteur latent ({M2} introns)")

    lp = np.empty((M2, S), dtype=np.float32)
    psi = np.empty((M2, S), dtype=np.float32)
    mu_all = np.empty((M2, S), dtype=np.float32)
    Bc_ok = np.empty((M2, S), dtype=bool)
    Kkeep = np.empty((M2, S), dtype=np.float32)
    Nkeep = np.empty((M2, S), dtype=np.float32)
    block = max(1000, min(20000, 4_000_000 // max(S, 1)))
    for a in range(0, M2, block):
        b = min(a + block, M2)
        cc = cols[a:b]
        kb0 = np.asarray(Dm[:, cc], dtype=np.float64).T
        bb0 = np.asarray(Bm[:, cc], dtype=np.float64).T
        nb0 = kb0 + bb0
        kb = np.rint(kb0)
        nb = np.maximum(np.rint(nb0), kb + 1)
        val = bb0 >= min_border          # couverture exploitable (tous)
        fit = val & ref[None, :]         # ... et sert de reference
        mu = robust_expected(kb, nb, X, valid=fit)
        msk = trim_mask(kb, nb, mu, valid=fit)
        rho = np.maximum(refine_rho(kb, nb, mu,
                                    rho_moments(kb, nb, mu, mask=msk),
                                    mask=msk), rho_floor)
        l, _ = bb_logp_screened(kb, nb, mu, np.repeat(rho[:, None], S, axis=1))
        lp[a:b] = l.reshape(b - a, S).astype(np.float32)
        mu_all[a:b] = mu.astype(np.float32)
        psi[a:b] = (kb0 / np.maximum(nb0, 1e-9)).astype(np.float32)
        Bc_ok[a:b] = val
        Kkeep[a:b] = kb0.astype(np.float32)
        Nkeep[a:b] = nb0.astype(np.float32)

    delta = psi - mu_all
    with np.errstate(invalid="ignore"):
        # mediane de cohorte = celle de la REFERENCE
        med = np.nanmedian(np.where(Bc_ok & ref[None, :], psi, np.nan), axis=1)
    med = np.where(np.isfinite(med), med, 0.0)
    testable = Bc_ok
    n_cov = (testable & ref[None, :]).sum(axis=1)
    Kc, Nc = Kkeep, Nkeep

    return dict(d=d, cols=cols, lp=lp, psi=psi, mu=mu_all, delta=delta,
                testable=testable, K=Kkeep, N=Nkeep, med=med, n_cov=n_cov,
                ref=ref, samples=samples)


def coverage_hits(wd, ann, padj_cutoff: float = 0.05,
                  delta_cutoff: float = 0.1, min_border: float = 10.0,
                  min_ratio: float = 0.1, rho_floor: float = 0.005,
                  min_samples: int = 25, top_n: int = 300,
                  top_per_sample: int = 20,
                  n_factors: Optional[int] = None,
                  samples_keep: Optional[Sequence[str]] = None,
                  model: Optional[dict] = None) -> List[dict]:
    """Evenements de retention vus par la couverture (voir coverage_model)."""
    from .stats import bh

    if model is None:
        model = coverage_model(wd, min_border=min_border, rho_floor=rho_floor,
                               min_samples=min_samples, n_factors=n_factors)
    if model is None:
        return []
    d, cols, lp, psi = model["d"], model["cols"], model["lp"], model["psi"]
    mu_all, delta, testable = model["mu"], model["delta"], model["testable"]
    Kc, Nc, med, n_cov = model["K"], model["N"], model["med"], model["n_cov"]
    ref, samples = model["ref"], model["samples"]
    M2 = len(cols)

    gname = list(ann.gene_name)
    gid = list(ann.gene_id)
    strand = ann.gene_strand
    # nombre d'introns testes PAR GENE : sert a corriger a l'interieur du gene
    gcol = np.asarray(d["gene"])[cols]
    n_by_gene = np.bincount(np.maximum(gcol, 0))
    n_g = n_by_gene[np.maximum(gcol, 0)].astype(np.float64)
    hits: List[dict] = []
    skipped = 0
    for si, sid in enumerate(samples):
        if samples_keep and sid not in samples_keep:
            continue
        p = np.exp(lp[:, si].astype(np.float64))
        p[~testable[:, si]] = np.nan
        padj = bh(p)
        # correction A L'INTERIEUR DU GENE (Bonferroni sur ses introns) :
        # c'est la lecture pertinente quand on regarde un gene candidat, mais
        # elle ne controle PAS le risque a l'echelle du genome.
        with np.errstate(invalid="ignore"):
            pg = np.minimum(p * n_g, 1.0)
        base = ((delta[:, si] >= delta_cutoff) & (psi[:, si] >= min_ratio)
                & testable[:, si])
        keep = base & (padj <= padj_cutoff)
        # liste de relecture : les N plus significatifs du sujet, meme si le
        # FDR genome entier ne les retient pas (un evenement a p = 5e-5 sur
        # 60 000 introns n'atteint jamais padj <= 0.05, et c'est pourtant le
        # genre de signal qu'un diagnostic doit voir passer)
        if top_per_sample > 0:
            cand = np.where(base & np.isfinite(p))[0]
            if len(cand):
                cand = cand[np.argsort(p[cand])[:top_per_sample]]
                keep[cand] = True
        sel = np.where(keep)[0]
        if len(sel) > top_n:
            skipped += len(sel) - top_n
            sel = sel[np.argsort(lp[sel, si])[:top_n]]
        for r in sel:
            c = int(cols[r])
            g = int(d["gene"][c])
            hits.append(dict(
                sample=sid, gene=gname[g] if 0 <= g < len(gname) else "",
                gene_id=gid[g] if 0 <= g < len(gid) else "",
                chrom=str(d["chrom"][c]),
                start=int(d["start"][c]) + 1, end=int(d["end"][c]),
                strand={1: "+", -1: "-"}.get(int(strand[g]) if 0 <= g < len(strand)
                                             else 0, "."),
                type="intron_retention", metrics="coverage", metric="coverage",
                k=int(round(Kc[r, si])), n=int(round(Nc[r, si])),
                psi=round(float(psi[r, si]), 4),
                mu=round(float(mu_all[r, si]), 4),
                delta=round(float(delta[r, si]), 4),
                psi_median_cohorte=round(float(med[r]), 4),
                p=float(p[r]), padj=float(padj[r]),
                mlog10p=float(-lp[r, si] / np.log(10)),
                padj_gene=float(pg[r]),
                retenu_par=("FDR" if padj[r] <= padj_cutoff else
                            ("gene" if pg[r] <= padj_cutoff else "candidate")),
                intron_len=int(d["length"][c]),
                reference=("yes" if ref[si] else "no"),
                n_samples_couverts=int(n_cov[r]), cov_row=int(r)))
    log(f"couverture : {len(hits)} evenement(s) sur {M2} introns testes"
        + (f" ({skipped} au-dela de --top)" if skipped else ""))
    return hits


def gene_coverage_table(wd, ann, gene: str,
                        min_border: float = 10.0) -> List[dict]:
    """Profondeurs brutes intron/bordures d'un gene, tous echantillons, sans seuil."""
    d = load_coverage(wd)
    if d is None:
        return []
    want = gene.strip().upper()
    gsel = {i for i, g in enumerate(ann.gene_name) if str(g).upper() == want}
    if not gsel:
        warn(f"gene inconnu dans l'annotation : {gene}")
        return []
    cols = [c for c in range(d["D"].shape[1]) if int(d["gene"][c]) in gsel]
    if not cols:
        warn(f"{gene} : aucun intron mesure")
        return []
    samples = [str(s) for s in d["samples"]]
    from .store import read_reference
    ref = read_reference(wd, samples)
    if ref.sum() < 5:
        ref = np.ones(len(samples), dtype=bool)

    # test beta-binomial restreint AUX INTRONS DE CE GENE : quand on sait quel
    # gene regarder, la correction multiple ne porte que sur ses introns, pas
    # sur les 60 000 du genome. C'est ce qui rend un evenement documentable.
    from .stats import (bb_logp_screened, bh, design_matrix, refine_rho,
                        rho_moments, robust_expected, trim_mask)
    Kg = np.asarray(d["D"][:, cols], dtype=np.float64).T
    Bg = np.asarray(d["Bw"][:, cols], dtype=np.float64).T
    Ng = Kg + Bg
    val = Bg >= min_border
    fit = val & ref[None, :]
    X = design_matrix(np.zeros((len(samples), 0)))
    kk = np.rint(Kg)
    nn = np.maximum(np.rint(Ng), kk + 1)
    mu_g = robust_expected(kk, nn, X, valid=fit)
    msk = trim_mask(kk, nn, mu_g, valid=fit)
    rho_g = np.maximum(refine_rho(kk, nn, mu_g,
                                  rho_moments(kk, nn, mu_g, mask=msk),
                                  mask=msk), 0.005)
    lpg, _ = bb_logp_screened(kk, nn, mu_g,
                              np.repeat(rho_g[:, None], len(samples), axis=1))
    lpg = lpg.reshape(Kg.shape)
    pg = np.exp(lpg)
    pg[~val] = np.nan
    padj_g = np.column_stack([bh(pg[:, i]) for i in range(len(samples))])

    out = []
    for ci, c in enumerate(cols):
        D = d["D"][:, c].astype(float)
        B1 = d["Bd"][:, c].astype(float)          # base bordante unique
        B = d["Bw"][:, c].astype(float)           # fenetre exonique bordante
        ok = B >= min_border                      # echantillons exploitables
        r = np.where(ok, D / np.maximum(B, 1e-9), np.nan)
        okr = ok & ref                            # ... et de reference
        if okr.any():
            med = float(np.nanmedian(r[okr]))
            sd = 1.4826 * float(np.nanmedian(np.abs(r[okr] - med))) + 0.01
        else:
            med, sd = float("nan"), 1.0
        # rang etabli UNIQUEMENT parmi les echantillons exploitables : un
        # rapport calcule sur 1 read de bordure n'est pas un resultat.
        for si, sid in enumerate(samples):
            rr = float(r[si]) if ok[si] else float("nan")
            out.append(dict(
                gene=gene, chrom=str(d["chrom"][c]),
                start=int(d["start"][c]) + 1, end=int(d["end"][c]),
                longueur=int(d["length"][c]), sample=sid,
                prof_intron=round(float(D[si]), 2),
                prof_bordure_1nt=round(float(B1[si]), 2),
                prof_bordures=round(float(B[si]), 2),
                exploitable="oui" if ok[si] else "non",
                reference="oui" if ref[si] else "non",
                ratio=(round(rr, 4) if ok[si] else ""),
                ratio_median_cohorte=round(med, 4) if med == med else "",
                z_robuste=(round((rr - med) / sd, 2) if ok[si] else ""),
                n_exploitables=int(okr.sum()),
                rang=(int(np.nansum(np.where(okr, r, np.nan) > rr)) + 1
                      if ok[si] else 0),
                psi=round(float(Kg[ci, si] / max(Ng[ci, si], 1e-9)), 4),
                psi_attendu=round(float(mu_g[ci, si]), 4),
                p_value=(float(pg[ci, si]) if ok[si] else ""),
                padj_gene=(float(padj_g[ci, si]) if ok[si] else "")))
    return out
