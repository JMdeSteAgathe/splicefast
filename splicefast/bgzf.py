"""
Lecture BGZF/BAM en Python pur : moteur de repli (--engine python) quand le
lecteur C++ sfscan n'est pas compile. Sert aussi a lire les en-tetes et index.

Contenu :
  * BgzfBlockStream : decompression sequentielle par blocs, avec offsets virtuels
  * BamHeader       : lecture de l'entete BAM (noms + longueurs des references)
  * read_bai        : index BAI -> points de decoupage pour la parallelisation
  * split_bam       : liste de tranches (voffset_debut, voffset_fin) equilibrees
"""
from __future__ import annotations

import os
import struct
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, List, Optional, Tuple

from .util import RobustReader, retry_io, winpath

BGZF_EOF = (
    b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43"
    b"\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
)


def voffset(coffset: int, uoffset: int) -> int:
    return (coffset << 16) | uoffset


def vsplit(v: int) -> Tuple[int, int]:
    return (v >> 16), (v & 0xFFFF)


class BgzfBlockStream:
    """Itere les blocs BGZF decompresses a partir d'un offset compresse donne.

    Renvoie des tuples (coffset_du_bloc, payload_decompresse).
    Ne charge jamais plus d'un bloc (64 Ko max) en memoire.
    """

    def __init__(self, path: str | os.PathLike, start_coffset: int = 0,
                 read_buffer: int = 4 << 20):
        self.path = winpath(path)
        self.fh = RobustReader(self.path, buffering=read_buffer)
        self.fh.seek(start_coffset)
        self.coffset = start_coffset
        self._decomp = zlib.decompressobj(-15)

    def close(self) -> None:
        try:
            self.fh.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()

    def blocks(self, stop_coffset: Optional[int] = None) -> Iterator[Tuple[int, bytes]]:
        fh = self.fh
        while True:
            if stop_coffset is not None and self.coffset > stop_coffset:
                return
            here = self.coffset
            head = fh.read(12)
            if len(head) < 12:
                return
            if head[0:2] != b"\x1f\x8b":
                raise ValueError(f"{self.path}: bloc BGZF invalide a l'offset {here}")
            xlen = struct.unpack_from("<H", head, 10)[0]
            extra = fh.read(xlen)
            if len(extra) < xlen:
                return
            bsize = None
            i = 0
            while i + 4 <= xlen:
                si1, si2, slen = struct.unpack_from("<BBH", extra, i)
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack_from("<H", extra, i + 4)[0] + 1
                    break
                i += 4 + slen
            if bsize is None:
                raise ValueError(f"{self.path}: champ BC absent (fichier non BGZF ?)")
            cdata_len = bsize - xlen - 20
            if not 0 <= cdata_len <= (1 << 18):
                raise ValueError(
                    f"{self.path}: en-tete BGZF incoherent a l'offset {here} "
                    f"(bsize={bsize}, xlen={xlen}) - fichier corrompu ?")
            cdata = fh.read(cdata_len)
            trailer = fh.read(8)
            if len(cdata) < cdata_len or len(trailer) < 8:
                return
            self.coffset = here + bsize
            if cdata_len == 0:
                # bloc EOF vide
                continue
            payload = zlib.decompress(cdata, -15)
            yield here, payload


@dataclass
class BamHeader:
    text: str
    ref_names: List[str]
    ref_lengths: List[int]
    data_start: int  # offset virtuel du premier enregistrement

    @property
    def n_ref(self) -> int:
        return len(self.ref_names)


def read_bam_header(path: str | os.PathLike) -> BamHeader:
    """Lit l'entete d'un BAM (peut chevaucher plusieurs blocs BGZF)."""
    buf = bytearray()
    blocks: List[Tuple[int, int]] = []  # (coffset, taille du payload)
    with BgzfBlockStream(path) as st:
        it = st.blocks()
        for coff, payload in it:
            blocks.append((coff, len(payload)))
            buf += payload
            if len(buf) < 12:
                continue
            if bytes(buf[:4]) != b"BAM\x01":
                raise ValueError(f"{path}: ce n'est pas un fichier BAM")
            l_text = struct.unpack_from("<i", buf, 4)[0]
            if len(buf) < 8 + l_text + 4:
                continue
            text = buf[8:8 + l_text].decode("utf-8", "replace")
            off = 8 + l_text
            n_ref = struct.unpack_from("<i", buf, off)[0]
            off += 4
            names: List[str] = []
            lengths: List[int] = []
            ok = True
            for _ in range(n_ref):
                if len(buf) < off + 4:
                    ok = False
                    break
                l_name = struct.unpack_from("<i", buf, off)[0]
                off += 4
                if len(buf) < off + l_name + 4:
                    ok = False
                    break
                names.append(buf[off:off + l_name - 1].decode("utf-8", "replace"))
                off += l_name
                lengths.append(struct.unpack_from("<i", buf, off)[0])
                off += 4
            if not ok:
                continue
            # convertir l'offset lineaire "off" en offset virtuel
            while True:
                rem = off
                for coff, size in blocks:
                    if rem < size:
                        return BamHeader(text, names, lengths, voffset(coff, rem))
                    rem -= size
                # l'entete finit pile en fin de bloc : lire le bloc suivant
                nxt = next(it, None)
                if nxt is None:
                    return BamHeader(text, names, lengths,
                                     voffset(st.coffset, 0))
                blocks.append((nxt[0], len(nxt[1])))
    raise ValueError(f"{path}: entete BAM tronquee")


def bai_path(bam: str | os.PathLike) -> Optional[Path]:
    p = Path(str(bam))
    for cand in (Path(str(p) + ".bai"), p.with_suffix(".bai"),
                 Path(str(p) + ".BAI"), p.with_suffix(".BAI")):
        if cand.exists():
            return cand
    return None


def read_bai_linear_offsets(path: str | os.PathLike) -> List[int]:
    """Renvoie tous les offsets virtuels de l'index lineaire (tries, uniques).

    Dans un BAM trie par coordonnees, ces offsets sont des debuts
    d'enregistrement garantis : ils constituent des points de decoupe surs.
    """
    def _slurp():
        with open(winpath(path), "rb") as fh:
            return fh.read()

    data = retry_io(_slurp, what=f"lecture de l'index {os.path.basename(str(path))}",
                    path=winpath(path))
    if data[:4] == b"BAI\x01":
        raw = data
    else:  # BAI parfois compresse en BGZF
        try:
            raw = zlib.decompress(data, 47)
        except Exception:
            raise ValueError(f"{path}: index BAI illisible")
        if raw[:4] != b"BAI\x01":
            raise ValueError(f"{path}: magic BAI absent")
    off = 4
    (n_ref,) = struct.unpack_from("<i", raw, off)
    off += 4
    offsets: List[int] = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", raw, off)
        off += 4
        for _ in range(n_bin):
            off += 4  # bin id
            (n_chunk,) = struct.unpack_from("<i", raw, off)
            off += 4
            off += 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", raw, off)
        off += 4
        if n_intv:
            arr = struct.unpack_from("<%dQ" % n_intv, raw, off)
            off += 8 * n_intv
            offsets.extend(v for v in arr if v)
    return sorted(set(offsets))


def read_bai_linear_index(path: str | os.PathLike) -> List[List[int]]:
    """Index lineaire par reference : lin[refid][pos >> 14] = offset virtuel.

    Permet d'aller lire directement une region donnee, sans parcourir tout le
    BAM. Les trous (0) sont combles par la valeur precedente, de sorte que
    lin[r][w] est toujours un point de depart sur (jamais apres la region).
    """
    def _slurp():
        with open(winpath(path), "rb") as fh:
            return fh.read()

    data = retry_io(_slurp, what=f"lecture de l'index {os.path.basename(str(path))}",
                    path=winpath(path))
    if data[:4] == b"BAI\x01":
        raw = data
    else:
        raw = zlib.decompress(data, 47)
        if raw[:4] != b"BAI\x01":
            raise ValueError(f"{path}: magic BAI absent")
    off = 4
    (n_ref,) = struct.unpack_from("<i", raw, off)
    off += 4
    out: List[List[int]] = []
    for _ in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", raw, off)
        off += 4
        for _ in range(n_bin):
            off += 4
            (n_chunk,) = struct.unpack_from("<i", raw, off)
            off += 4
            off += 16 * n_chunk
        (n_intv,) = struct.unpack_from("<i", raw, off)
        off += 4
        if n_intv:
            arr = list(struct.unpack_from("<%dQ" % n_intv, raw, off))
            off += 8 * n_intv
            last = 0
            for i, v in enumerate(arr):
                if v:
                    last = v
                else:
                    arr[i] = last
            out.append(arr)
        else:
            out.append([])
    return out


def split_bam(bam: str | os.PathLike, n_parts: int,
              min_chunk_bytes: int = 32 << 20) -> List[Tuple[int, int]]:
    """Decoupe le BAM en tranches [v_debut, v_fin) sur des frontieres
    d'enregistrement garanties (index lineaire du BAI).

    Renvoie [(v0, v1), (v1, v2), ...]; la derniere tranche va jusqu'a la fin
    (v_fin = 0 signifie 'jusqu'au bout du fichier').
    Si aucun BAI n'est disponible : une seule tranche (lecture sequentielle).
    """
    hdr = read_bam_header(bam)
    size = os.path.getsize(winpath(bam))
    if n_parts <= 1:
        return [(hdr.data_start, 0)]
    bi = bai_path(bam)
    if bi is None:
        return [(hdr.data_start, 0)]
    try:
        lin = read_bai_linear_offsets(bi)
    except Exception:
        return [(hdr.data_start, 0)]
    lin = [v for v in lin if v > hdr.data_start]
    if not lin:
        return [(hdr.data_start, 0)]
    # cible : n_parts tranches, mais au moins min_chunk_bytes de donnees compressees
    target = max(min_chunk_bytes, size // max(1, n_parts))
    cuts = [hdr.data_start]
    last_c = vsplit(hdr.data_start)[0]
    for v in lin:
        c = vsplit(v)[0]
        if c - last_c >= target:
            cuts.append(v)
            last_c = c
    out = [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]
    out.append((cuts[-1], 0))
    return out
