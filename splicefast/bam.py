"""
Extraction des comptages a partir d'un BAM, en Python pur.

Pour chaque tranche de BAM on produit :
  * les jonctions (operations 'N' du CIGAR)  -> reads coupes    (K)
  * les reads NON coupes chevauchant chaque site d'epissage    (retention)

Astuce de performance pour les reads non coupes : au lieu de tester chaque
site pour chaque read, on utilise un tableau de differences + searchsorted
vectorise -> cout O(1) par bloc aligne, quel que soit le nombre de sites.

Coordonnees : tout est en coordonnees GLOBALES 0-based (offset cumule des
references) ce qui permet de n'avoir qu'un seul tableau trie de sites et une
seule cle entiere par jonction.
"""
from __future__ import annotations

import struct
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .bgzf import BgzfBlockStream, read_bam_header, split_bam, vsplit

# flags exclus par defaut : unmapped(4) secondary(256) qcfail(512) dup(1024)
# supplementary(2048)
DEFAULT_EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400 | 0x800

_CORE = struct.Struct("<iiBBHHHi")
_U32 = struct.Struct("<I")

LEN_BITS = 24
LEN_MASK = (1 << LEN_BITS) - 1


def pack_junction(gstart: int, length: int) -> int:
    return (gstart << LEN_BITS) | (length & LEN_MASK)


def unpack_junctions(keys: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    gstart = keys >> LEN_BITS
    length = keys & LEN_MASK
    return gstart, length


def ref_offsets(ref_lengths: Sequence[int], pad: int = 1000) -> np.ndarray:
    """Offsets cumules par reference (coordonnees globales)."""
    off = np.zeros(len(ref_lengths) + 1, dtype=np.int64)
    acc = 0
    for i, L in enumerate(ref_lengths):
        off[i] = acc
        acc += int(L) + pad
    off[len(ref_lengths)] = acc
    return off


class ChunkResult:
    __slots__ = ("junc", "nonsplit", "n_reads", "n_spliced", "n_blocks")

    def __init__(self, junc, nonsplit, n_reads, n_spliced, n_blocks):
        self.junc: Dict[int, int] = junc
        self.nonsplit: Optional[np.ndarray] = nonsplit
        self.n_reads = n_reads
        self.n_spliced = n_spliced
        self.n_blocks = n_blocks


def scan_chunk(
    bam: str,
    vstart: int,
    vend: int,
    ref_off: np.ndarray,
    sites: Optional[np.ndarray] = None,
    anchor: int = 5,
    min_mapq: int = 0,
    exclude_flags: int = DEFAULT_EXCLUDE,
    min_intron: int = 20,
    max_intron: int = 1_000_000,
    min_anchor_split: int = 5,
    read_buffer: int = 4 << 20,
    flush_blocks: int = 400_000,
) -> ChunkResult:
    """Scanne [vstart, vend) d'un BAM et renvoie jonctions + reads non coupes.

    sites : positions GLOBALES triees des frontieres exon/intron (0-based ;
            la frontiere 'p' signifie entre la base p-1 et la base p).
    anchor: nombre de bases requises de chaque cote de la frontiere pour
            compter un read comme 'non coupe' (equivalent minAnchor de FRASER).
    """
    junc: Dict[int, int] = {}
    jget = junc.get
    n_reads = 0
    n_spliced = 0
    n_blocks_tot = 0

    have_sites = sites is not None and len(sites) > 0
    if have_sites:
        nsite = len(sites)
        diff = np.zeros(nsite + 2, dtype=np.int32)
        lo_lst: List[int] = []
        hi_lst: List[int] = []
        sites_lo = sites  # alias local
    else:
        nsite = 0
        diff = None

    start_c, start_u = vsplit(vstart)
    if vend:
        stop_c, stop_u = vsplit(vend)
    else:
        stop_c, stop_u = None, None

    buf = bytearray()
    stream = BgzfBlockStream(bam, start_c, read_buffer=read_buffer)
    first = True
    finished = False

    def flush_sites():
        nonlocal lo_lst, hi_lst
        if not lo_lst:
            return
        a = np.fromiter(lo_lst, dtype=np.int64, count=len(lo_lst))
        b = np.fromiter(hi_lst, dtype=np.int64, count=len(hi_lst))
        lo_lst = []
        hi_lst = []
        # frontieres s telles que  a <= s <= b
        i0 = np.searchsorted(sites_lo, a, side="left")
        i1 = np.searchsorted(sites_lo, b, side="right")
        keep = i1 > i0
        if not keep.any():
            return
        np.add.at(diff, i0[keep], 1)
        np.add.at(diff, i1[keep], -1)

    try:
        for coff, payload in stream.blocks(stop_coffset=stop_c):
            if first:
                first = False
                if stop_c is not None and coff == stop_c:
                    buf += payload[start_u:stop_u]
                    finished = True
                else:
                    buf += payload[start_u:]
            elif stop_c is not None and coff == stop_c:
                buf += payload[:stop_u]
                finished = True
            else:
                buf += payload

            # ---- boucle de parsing des enregistrements ----
            off = 0
            n = len(buf)
            while True:
                if n - off < 4:
                    break
                (bs,) = _U32.unpack_from(buf, off)
                if n - off < 4 + bs:
                    break
                (refid, pos, l_rn, mapq, _bin, n_cig, flag,
                 _l_seq) = _CORE.unpack_from(buf, off + 4)
                off_next = off + 4 + bs
                if refid < 0 or (flag & exclude_flags) or mapq < min_mapq:
                    off = off_next
                    continue
                n_reads += 1
                gpos = ref_off[refid] + pos

                if n_cig == 0:
                    off = off_next
                    continue
                if n_cig == 1:
                    (c,) = _U32.unpack_from(buf, off + 36 + l_rn)
                    if (c & 0xF) in (0, 7, 8):
                        ln = c >> 4
                        if have_sites and ln >= 2 * anchor:
                            lo_lst.append(gpos + anchor)
                            hi_lst.append(gpos + ln - anchor)
                            n_blocks_tot += 1
                    off = off_next
                    continue

                cig_off = off + 36 + l_rn
                cur = gpos
                blk_start = gpos
                spliced = False
                pend = 0  # jonction en attente de validation de son ancrage droit
                for i in range(n_cig):
                    (c,) = _U32.unpack_from(buf, cig_off + 4 * i)
                    op = c & 0xF
                    ln = c >> 4
                    if op == 3:  # N -> jonction
                        blk_len = cur - blk_start
                        # la jonction precedente voit ici son ancrage droit
                        if pend and blk_len >= min_anchor_split:
                            junc[pend] = (jget(pend) or 0) + 1
                            spliced = True
                        if blk_len >= min_anchor_split and \
                                min_intron <= ln <= max_intron:
                            pend = pack_junction(cur, ln)
                        else:
                            pend = 0
                        if have_sites and blk_len >= 2 * anchor:
                            lo_lst.append(blk_start + anchor)
                            hi_lst.append(cur - anchor)
                            n_blocks_tot += 1
                        cur += ln
                        blk_start = cur
                    elif op in (0, 7, 8, 2):  # M/=/X/D : consomme la reference
                        cur += ln
                    # I, S, H, P : ne consomment pas la reference
                if pend and (cur - blk_start) >= min_anchor_split:
                    junc[pend] = (jget(pend) or 0) + 1
                    spliced = True
                if have_sites and (cur - blk_start) >= 2 * anchor:
                    lo_lst.append(blk_start + anchor)
                    hi_lst.append(cur - anchor)
                    n_blocks_tot += 1
                if spliced:
                    n_spliced += 1
                if have_sites and len(lo_lst) >= flush_blocks:
                    flush_sites()
                off = off_next

            if off:
                del buf[:off]
            if finished:
                break
    finally:
        stream.close()

    nonsplit = None
    if have_sites:
        flush_sites()
        nonsplit = np.cumsum(diff[:nsite], dtype=np.int64).astype(np.int32)
    return ChunkResult(junc, nonsplit, n_reads, n_spliced, n_blocks_tot)


def scan_chunk_task(args) -> Tuple[Dict[int, int], Optional[np.ndarray], int, int]:
    """Point d'entree picklable pour ProcessPoolExecutor (Windows = spawn)."""
    (bam, vstart, vend, ref_off, sites, opts) = args
    res = scan_chunk(bam, vstart, vend, ref_off, sites, **opts)
    return res.junc, res.nonsplit, res.n_reads, res.n_spliced


def bam_chunks(bam: str, n_parts: int, min_chunk_bytes: int = 32 << 20):
    return split_bam(bam, n_parts, min_chunk_bytes=min_chunk_bytes)


def bam_refs(bam: str) -> Tuple[List[str], List[int]]:
    hdr = read_bam_header(bam)
    return hdr.ref_names, hdr.ref_lengths


def bam_is_paired(bam: str, n: int = 2000) -> bool:
    """Teste rapidement si la librairie est paired-end (1er bloc du BAM)."""
    hdr = read_bam_header(bam)
    c, u = vsplit(hdr.data_start)
    buf = bytearray()
    paired = 0
    seen = 0
    with BgzfBlockStream(bam, c) as st:
        for i, (coff, payload) in enumerate(st.blocks()):
            buf += payload[u:] if i == 0 else payload
            off = 0
            while len(buf) - off >= 4:
                (bs,) = _U32.unpack_from(buf, off)
                if len(buf) - off < 4 + bs:
                    break
                (_r, _p, _l, _m, _b, _nc, flag, _ls) = _CORE.unpack_from(buf, off + 4)
                seen += 1
                if flag & 0x1:
                    paired += 1
                off += 4 + bs
                if seen >= n:
                    return paired > seen * 0.5
            del buf[:off]
            if i > 40:
                break
    return paired > max(1, seen) * 0.5
