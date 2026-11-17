"""
Auto-test : fabrique un mini run RNA-seq simule (BAM + GTF), y injecte des
evenements aberrants connus, puis verifie que le pipeline les retrouve.

Utile pour valider l'installation en 1 minute, sans aucune donnee reelle :
    python3 -m splicefast selftest -w /tmp/splicefast_selftest
"""
from __future__ import annotations

import struct
import zlib
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .util import log, warn, winpath

# --------------------------------------------------------------------------- #
# ecriture BAM + BAI (uniquement pour l'auto-test)
# --------------------------------------------------------------------------- #

def _bgzf_block(data: bytes) -> bytes:
    co = zlib.compressobj(6, zlib.DEFLATED, -15)
    comp = co.compress(data) + co.flush()
    bsize = len(comp) + 25 + 1
    if bsize > 65536:
        raise ValueError("bloc trop grand")
    head = struct.pack("<BBBBIBBH", 31, 139, 8, 4, 0, 0, 255, 6)
    extra = struct.pack("<BBHH", 66, 67, 2, bsize - 1)
    tail = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return head + extra + comp + tail


EOF_BLOCK = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000")

SEQ_NT = {"A": 1, "C": 2, "G": 4, "T": 8, "N": 15}


def _pack_seq(seq: str) -> bytes:
    out = bytearray()
    for i in range(0, len(seq), 2):
        hi = SEQ_NT.get(seq[i], 15)
        lo = SEQ_NT.get(seq[i + 1], 15) if i + 1 < len(seq) else 0
        out.append((hi << 4) | lo)
    return bytes(out)


def _record(name: str, refid: int, pos: int, cigar: Sequence[Tuple[int, int]],
            seqlen: int, flag: int = 0, mapq: int = 60) -> bytes:
    nm = name.encode() + b"\x00"
    cig = b"".join(struct.pack("<I", (ln << 4) | op) for op, ln in cigar)
    seq = _pack_seq("A" * seqlen)
    qual = bytes([30] * seqlen)
    body = (struct.pack("<iiBBHHHiiii", refid, pos, len(nm), mapq, 0,
                        len(cigar), flag, seqlen, -1, -1, 0)
            + nm + cig + seq + qual)
    return struct.pack("<I", len(body)) + body


class BamWriter:
    """Ecrivain BAM minimal + index BAI (assez pour splicefast)."""

    def __init__(self, path: Path, refs: Sequence[Tuple[str, int]]):
        self.path = Path(path)
        self.refs = list(refs)
        self.fh = open(winpath(self.path), "wb")
        text = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
            f"@SQ\tSN:{n}\tLN:{L}\n" for n, L in refs)
        hdr = b"BAM\x01" + struct.pack("<i", len(text)) + text.encode()
        hdr += struct.pack("<i", len(refs))
        for n, L in refs:
            nb = n.encode() + b"\x00"
            hdr += struct.pack("<i", len(nb)) + nb + struct.pack("<i", L)
        for i in range(0, len(hdr), 60000):
            self.fh.write(_bgzf_block(hdr[i:i + 60000]))
        self.buf = bytearray()
        self.coffset = self.fh.tell()
        # index lineaire : voffset du 1er read de chaque fenetre de 16 kb
        self.lin: Dict[int, Dict[int, int]] = {i: {} for i in range(len(refs))}
        self.first: Dict[int, int] = {}
        self.last: Dict[int, int] = {}

    def _flush(self, force: bool = False) -> None:
        while len(self.buf) >= 60000 or (force and self.buf):
            take = min(60000, len(self.buf))
            self.fh.write(_bgzf_block(bytes(self.buf[:take])))
            del self.buf[:take]
            self.coffset = self.fh.tell()

    def write(self, refid: int, pos: int, rec: bytes) -> None:
        v = (self.coffset << 16) | len(self.buf)
        w = pos >> 14
        self.lin[refid].setdefault(w, v)
        self.first.setdefault(refid, v)
        self.last[refid] = v
        self.buf += rec
        self._flush()

    def close(self) -> None:
        self._flush(force=True)
        self.fh.write(EOF_BLOCK)
        self.fh.close()
        self._write_bai()

    def _write_bai(self) -> None:
        out = bytearray(b"BAI\x01" + struct.pack("<i", len(self.refs)))
        for i, (n, L) in enumerate(self.refs):
            lin = self.lin[i]
            if not lin:
                out += struct.pack("<i", 0) + struct.pack("<i", 0)
                continue
            # un seul bin (0 = toute la reference) avec un chunk global
            out += struct.pack("<i", 1)
            out += struct.pack("<I", 0) + struct.pack("<i", 1)
            out += struct.pack("<QQ", self.first[i], (self.last[i] >> 16 << 16)
                               + 65535)
            nintv = max(lin) + 1
            out += struct.pack("<i", nintv)
            cur = self.first[i]
            for w in range(nintv):
                cur = lin.get(w, cur)
                out += struct.pack("<Q", cur)
        with open(winpath(Path(str(self.path) + ".bai")), "wb") as fh:
            fh.write(bytes(out))


# --------------------------------------------------------------------------- #
# jeu de donnees simule
# --------------------------------------------------------------------------- #

def simulate(root: Path, n_samples: int = 12, n_genes: int = 24,
             depth: int = 45, read_len: int = 100, seed: int = 3):
    rng = np.random.default_rng(seed)
    root.mkdir(parents=True, exist_ok=True)
    bamdir = root / "bam"
    bamdir.mkdir(exist_ok=True)
    chrom = "chr1"
    exon_len = 180
    genes = []
    pos = 5000
    for g in range(n_genes):
        exons = []
        for e in range(5):
            exons.append((pos, pos + exon_len))
            pos += exon_len + int(rng.integers(600, 2500))
        genes.append(exons)
        pos += 3000
    ref_len = pos + 10000

    # GTF
    gtf = root / "annotation.gtf"
    with open(winpath(gtf), "w", encoding="utf-8") as fh:
        for g, exons in enumerate(genes):
            gid = f"GENE{g:03d}"
            fh.write(f"{chrom}\tsim\tgene\t{exons[0][0]+1}\t{exons[-1][1]}\t.\t+\t.\t"
                     f'gene_id "{gid}"; gene_name "{gid}"; gene_type "protein_coding";\n')
            for e, (s, en) in enumerate(exons):
                fh.write(f"{chrom}\tsim\texon\t{s+1}\t{en}\t.\t+\t.\t"
                         f'gene_id "{gid}"; transcript_id "{gid}.1"; '
                         f'gene_name "{gid}"; gene_type "protein_coding"; '
                         f'exon_number {e+1};'
                         + (' tag "Ensembl_canonical";' if g % 2 == 0 else '') + '\n')

    # evenements injectes
    spikes = {
        3: ("exon_skipping", 5),      # echantillon 3, gene 5 : saut d'exon 3
        5: ("intron_retention", 8),   # echantillon 5, gene 8 : retention intron 2
        7: ("pseudo_exon", 12),       # echantillon 7, gene 12 : pseudo-exon
        9: ("novel_donor", 17),       # echantillon 9, gene 17 : donneur nouveau
    }
    expected = []
    samples = []
    for s in range(n_samples):
        sid = f"SIM{s:02d}"
        bam = bamdir / f"{sid}.bam"
        w = BamWriter(bam, [(chrom, ref_len)])
        recs: List[Tuple[int, bytes]] = []
        rid = 0
        # effet "batch" : deux groupes de profondeur differente
        scale = 1.0 if s % 2 == 0 else 0.6
        for g, exons in enumerate(genes):
            for j in range(len(exons) - 1):
                istart, iend = exons[j][1], exons[j + 1][0]
                nreads = max(4, int(rng.poisson(depth * scale)))
                skip_this = (s in spikes and spikes[s][1] == g
                             and spikes[s][0] == "exon_skipping" and j == 1)
                ir_this = (s in spikes and spikes[s][1] == g
                           and spikes[s][0] == "intron_retention" and j == 1)
                pex_this = (s in spikes and spikes[s][1] == g
                            and spikes[s][0] == "pseudo_exon" and j == 1)
                nd_this = (s in spikes and spikes[s][1] == g
                           and spikes[s][0] == "novel_donor" and j == 1)
                if ir_this:
                    nreads = max(3, nreads // 5)
                if pex_this or nd_this:
                    nreads = max(3, nreads // 3)
                # reads epissés normaux
                for _ in range(nreads):
                    off = int(rng.integers(10, read_len - 10))
                    p = istart - off
                    recs.append((p, _record(f"r{rid}", 0, p,
                                            [(0, off), (3, iend - istart),
                                             (0, read_len - off)], read_len)))
                    rid += 1
                # reads couvrant le CORPS de l'intron (pre-ARNm) : fond
                # faible partout, massif en cas de retention -> c'est ce que
                # mesure le module de couverture
                n_in = int(rng.poisson(2 * scale)) if not ir_this \
                    else int(rng.poisson(45))
                for _ in range(n_in):
                    hi = max(istart + 1, iend - read_len)
                    p_in = int(rng.integers(istart, hi))
                    recs.append((p_in, _record(f"i{rid}", 0, p_in,
                                               [(0, read_len)], read_len)))
                    rid += 1
                # reads non coupés (fond : bord de l'intron)
                bg = max(1, int(rng.poisson(3 * scale)))
                if ir_this:
                    bg = int(rng.poisson(50))
                for _ in range(bg):
                    for b in (istart, iend):
                        off = int(rng.integers(20, read_len - 20))
                        p = b - off
                        recs.append((p, _record(f"n{rid}", 0, p,
                                                [(0, read_len)], read_len)))
                        rid += 1
                # evenements aberrants
                if skip_this:
                    tgt = exons[j + 2][0]
                    for _ in range(int(depth * 0.8)):
                        off = int(rng.integers(10, read_len - 10))
                        p = istart - off
                        recs.append((p, _record(f"s{rid}", 0, p,
                                                [(0, off), (3, tgt - istart),
                                                 (0, read_len - off)], read_len)))
                        rid += 1
                    expected.append((sid, "exon_skipping", chrom, istart + 1, tgt))
                if pex_this:
                    a = istart + 400
                    d = a + 120
                    for _ in range(int(depth * 0.7)):
                        off = int(rng.integers(10, read_len - 10))
                        p = istart - off
                        recs.append((p, _record(f"p{rid}", 0, p,
                                                [(0, off), (3, a - istart),
                                                 (0, read_len - off)], read_len)))
                        rid += 1
                        off = int(rng.integers(10, read_len - 10))
                        p = d - off
                        recs.append((p, _record(f"q{rid}", 0, p,
                                                [(0, off), (3, iend - d),
                                                 (0, read_len - off)], read_len)))
                        rid += 1
                    expected.append((sid, "pseudo_exon", chrom, istart + 1, a))
                if nd_this:
                    nd = istart + 60          # decalage bien au-dela de 20 pb
                    for _ in range(int(depth * 0.8)):
                        off = int(rng.integers(10, read_len - 10))
                        p = nd - off
                        recs.append((p, _record(f"d{rid}", 0, p,
                                                [(0, off), (3, iend - nd),
                                                 (0, read_len - off)], read_len)))
                        rid += 1
                    expected.append((sid, "novel_donor", chrom, nd + 1, iend))
                if ir_this:
                    expected.append((sid, "intron_retention", chrom, istart + 1,
                                     istart + 1))
        recs.sort(key=lambda x: x[0])
        for p, r in recs:
            w.write(0, p, r)
        w.close()
        samples.append((sid, str(bam)))
    return dict(gtf=str(gtf), samples=samples, expected=expected,
                bamdir=str(bamdir))


def compare_engines(root: Path, sim: dict, ann, ref_names, ref_lengths,
                    ref_off, sites_path: str, h: str, opts: dict,
                    main_wd) -> bool:
    """Le lecteur C++ doit donner EXACTEMENT les memes nombres que le Python."""
    from .count import count_all
    from .coverage import build_cov_targets, measure_coverage
    from .scanner import find_sfscan
    from .util import Workdir

    if find_sfscan() is None:
        warn("sfscan non compile : comparaison C++ / Python sautee")
        return True
    sub = sim["samples"][:4]
    wd2 = Workdir(root / "work_python").create()
    count_all(wd2, sub, ref_off, sites_path, h, opts, threads=1, force=True,
              engine="python")
    tg = build_cov_targets(ann, ref_names, ref_lengths)
    measure_coverage(wd2, sub, tg, ref_lengths, threads=1, engine="python")
    ok = True
    for sid, _ in sub:
        a = np.load(winpath(main_wd.sample_counts(sid)), allow_pickle=True)
        b = np.load(winpath(wd2.sample_counts(sid)), allow_pickle=True)
        for k in ("jkey", "jcount", "nonsplit"):
            ok &= bool(np.array_equal(a[k], b[k]))
        ok &= bool(np.array_equal(a["meta"][:3], b["meta"][:3]))
        ca = np.load(winpath(main_wd.covcounts / f"{sid}.npz"))["res"]
        cb = np.load(winpath(wd2.covcounts / f"{sid}.npz"))["res"]
        ok &= bool(np.array_equal(ca, cb))
    log(f"  lecteur C++ vs Python sur {len(sub)} echantillons : "
        f"{'IDENTIQUES' if ok else 'DIFFERENTS !'}")
    return ok


def run_selftest(workdir: Path, threads: int = 2) -> bool:
    from .analyze import analyze
    from .cli import _project
    from .count import count_all
    from .gtf import Annotation
    from .report import make_reports
    from .store import build_store
    from .util import Workdir

    root = Path(workdir)
    log(f"Auto-test dans {root}")
    sim = simulate(root / "sim")
    wd = Workdir(root / "work").create()
    ann, ref_names, ref_lengths, ref_off, sites_path, h = _project(
        wd, sim["samples"], sim["gtf"], force_ann=True)
    opts = dict(anchor=5, min_mapq=0, min_intron=20, max_intron=1000000,
                min_anchor_split=5)
    count_all(wd, sim["samples"], ref_off, sites_path, h, opts, threads=threads,
              force=True)
    build_store(wd, sim["samples"], ann, ref_names, ref_lengths,
                min_total=5, min_samples=1, block_size=5000)
    analyze(wd, metrics=("jaccard", "psi5", "psi3", "theta"), n_factors=3,
            sub_rows=5000, min_k_one=10, min_n_quantile=10, min_cov=5)
    # module dedie de retention (couverture des introns)
    from .coverage import (build_cov_targets, measure_coverage, save_coverage)
    tg = build_cov_targets(ann, ref_names, ref_lengths)
    S_, M_ = measure_coverage(wd, sim["samples"], tg, ref_lengths,
                              threads=threads)
    save_coverage(wd, tg, sim["samples"], S_, M_)

    hits = make_reports(wd, ann, padj_cutoff=0.1, delta_cutoff=0.05, top_n=100)

    found = 0
    log("--- verification des evenements injectes ---")
    for sid, etype, chrom, start, end in sim["expected"]:
        ok = [h for h in hits if h["sample"] == sid and h["chrom"] == chrom
              and abs(h["start"] - start) <= 2
              and (etype == "intron_retention" or abs((h["end"] or 0) - end) <= 2)]
        types = {x.get("type") for x in ok}
        good = bool(ok) and (etype in types or (
            etype == "novel_donor" and types & {"novel_donor", "splice_site_shift"}))
        found += int(good)
        log(f"  {sid} {etype} {chrom}:{start}-{end} -> "
            f"{'TROUVE' if good else 'MANQUE'} {sorted(types) if ok else ''}")
    # le module de couverture doit voir la retention massive
    cov = [h for h in hits if h.get("metrics") == "coverage"]
    ir = [h for h in cov if h["sample"] == "SIM05"]
    log(f"  module couverture : {len(cov)} evenement(s), "
        f"retention SIM05 {'TROUVEE' if ir else 'MANQUEE'}"
        + (f" (ratio {ir[0]['psi']} vs {ir[0]['psi_median_cohorte']} "
           f"dans la cohorte)" if ir else ""))
    same = compare_engines(root, sim, ann, ref_names, ref_lengths, ref_off,
                           sites_path, h, opts, wd)
    log(f"Auto-test : {found}/{len(sim['expected'])} evenements retrouves, "
        f"{len(hits)} evenements rapportes au total")
    log(f"Visualiseur : {wd.results / 'SpliceFast_viewer.html'}")
    return found >= len(sim["expected"]) - 1 and bool(ir) and same
