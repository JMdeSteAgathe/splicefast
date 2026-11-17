"""
`splicefast probe` : relire une region DIRECTEMENT dans le BAM et comparer a ce
que SpliceFast a stocke.

    python -m splicefast probe -w DOSSIER --sample ECH02-PAX_S61 \
        --region chr5:37025000-37028000

Pour chaque jonction de la region, on affiche :
  * le nombre de reads dans le BAM, ventile par categorie de filtre
    (retenu / doublon PCR / secondaire / supplementaire / qcfail / MAPQ) ;
  * le nombre stocke dans le fichier de comptage de l'echantillon
    (counts/<echantillon>.npz) ;
  * la valeur dans la matrice de la cohorte (store/) et la mediane chez les
    autres echantillons.

Et, pour chaque site d'epissage annote de la region, les reads NON coupes
(ceux qui alimentent theta) dans le BAM et dans le fichier de comptage.

C'est la commande a lancer quand un evenement vu dans IGV / un viewer sashimi
/ FRASER n'apparait pas : elle dit en une fois si les reads sont filtres
(doublons, MAPQ...), si le fichier de comptage ne correspond pas au BAM
(comptage perime, mauvais fichier), ou si le probleme est en aval (modele,
seuils -> voir `why`).
"""
from __future__ import annotations

import os
import re
import struct
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .bam import DEFAULT_EXCLUDE, pack_junction, ref_offsets
from .bgzf import (BgzfBlockStream, bai_path, read_bai_linear_index,
                   read_bam_header, vsplit)
from .util import log, warn, winpath

_CORE = struct.Struct("<iiBBHHHi")
_U32 = struct.Struct("<I")

CLASSES = ("kept", "duplicate", "secondary", "supplementary", "qcfail",
           "low_mapq")


def parse_region(text: str) -> Tuple[str, int, int]:
    m = re.match(r"^\s*([^:\s]+)\s*:\s*([\d,]+)\s*-\s*([\d,]+)\s*$", text)
    if not m:
        raise SystemExit(f"region illisible : {text!r} (attendu chr5:37025000-37028000)")
    return m.group(1), int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))


def _klass(flag: int, mapq: int, min_mapq: int) -> Optional[str]:
    if flag & 0x4:
        return None
    if flag & 0x100:
        return "secondary"
    if flag & 0x800:
        return "supplementary"
    if flag & 0x200:
        return "qcfail"
    if flag & 0x400:
        return "duplicate"
    if mapq < min_mapq:
        return "low_mapq"
    return "kept"


def scan_region(bam: str, chrom: str, start1: int, end1: int,
                anchor: int = 5, min_anchor_split: int = 5,
                min_intron: int = 20, max_intron: int = 1_000_000,
                min_mapq: int = 0, sites_local: Optional[np.ndarray] = None):
    """Lit les reads chevauchant [start1, end1] (1-based) d'un BAM.

    Renvoie (junctions, nonsplit, n_reads) ou
      junctions[(intron_start0, length)][classe] = nombre de reads
      nonsplit[site0][classe] = reads non coupes au site (frontiere entre
                                 la base site0-1 et site0, 0-based)
    """
    hdr = read_bam_header(bam)
    names = list(hdr.ref_names)
    if chrom not in names:
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        if alt not in names:
            raise SystemExit(f"{chrom} absent de l'en-tete de {bam}")
        chrom = alt
    rid = names.index(chrom)
    s0, e0 = start1 - 1, end1
    bi = bai_path(bam)
    if bi is None:
        warn(f"pas d'index .bai pour {bam} : lecture depuis le debut (lent)")
        vstart = hdr.data_start
    else:
        lin = read_bai_linear_index(bi)[rid]
        if not lin:
            return {}, {}, 0
        # les reads longs (introns) peuvent commencer bien avant la region
        w = max(0, (s0 - 1_000_000) >> 14)
        w = min(w, len(lin) - 1)
        vstart = lin[w]
    coff, uoff = vsplit(vstart)
    junc: Dict[Tuple[int, int], Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    nons: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    sites = (np.asarray(sorted(set(int(x) for x in sites_local)), dtype=np.int64)
             if sites_local is not None and len(sites_local) else np.zeros(0, np.int64))
    n_reads = 0
    buf = bytearray()
    first = True
    done = False
    with BgzfBlockStream(bam, coff) as st:
        for _c, payload in st.blocks():
            buf += payload[uoff:] if first else payload
            first = False
            off, n = 0, len(buf)
            while n - off >= 4:
                (bs,) = _U32.unpack_from(buf, off)
                if n - off < 4 + bs:
                    break
                (refid, pos, l_rn, mapq, _b, n_cig, flag, _l) = _CORE.unpack_from(buf, off + 4)
                nxt = off + 4 + bs
                off = nxt
                if refid != rid:
                    if refid > rid or refid < 0:
                        done = True
                        break
                    continue
                if pos > e0:
                    done = True
                    break
                k = _klass(flag, mapq, min_mapq)
                if k is None or n_cig == 0:
                    continue
                cig = nxt - bs + 32 + l_rn
                cur = pos
                blk = pos
                blocks = []
                juncs = []
                pend = None
                for i in range(n_cig):
                    (c,) = _U32.unpack_from(buf, cig + 4 * i)
                    op, ln = c & 0xF, c >> 4
                    if op == 3:
                        bl = cur - blk
                        if pend is not None and bl >= min_anchor_split:
                            juncs.append(pend)
                        pend = ((cur, ln) if bl >= min_anchor_split
                                and min_intron <= ln <= max_intron else None)
                        blocks.append((blk, cur))
                        cur += ln
                        blk = cur
                    elif op in (0, 7, 8, 2):
                        cur += ln
                if pend is not None and cur - blk >= min_anchor_split:
                    juncs.append(pend)
                blocks.append((blk, cur))
                if cur < s0:
                    continue
                n_reads += 1
                for j in juncs:
                    if j[0] + j[1] >= s0 and j[0] <= e0:
                        junc[j][k] += 1
                if len(sites):
                    for a, b in blocks:
                        if b - a < 2 * anchor:
                            continue
                        lo, hi = a + anchor, b - anchor
                        i0 = np.searchsorted(sites, lo, side="left")
                        i1 = np.searchsorted(sites, hi, side="right")
                        for s in sites[i0:i1].tolist():
                            nons[s][k] += 1
            del buf[:off]
            if done:
                break
    return junc, nons, n_reads


def probe(wd, ann, sample: Optional[str], region: str, bam: Optional[str] = None,
          min_mapq: int = 0, anchor: int = 5) -> List[dict]:
    from .store import Store
    from .util import safe_name

    chrom, a1, b1 = parse_region(region)
    st = None
    try:
        st = Store(wd)
    except Exception:  # noqa: BLE001
        warn("store absent : comparaison au BAM et au fichier de comptage seulement")
    if bam is None:
        if st is None or sample not in st.samples:
            raise SystemExit(f"echantillon inconnu : {sample} (ou passer --bam)")
        bam = st.meta["bams"][st.samples.index(sample)]
    log(f"BAM lu : {bam}")
    if not os.path.exists(winpath(bam)):
        raise SystemExit(f"BAM introuvable : {bam}")

    ref_names = st.ref_names if st is not None else list(read_bam_header(bam).ref_names)
    ref_lengths = (st.ref_lengths if st is not None
                   else list(read_bam_header(bam).ref_lengths))
    if chrom not in ref_names:
        alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
        chrom = alt if alt in ref_names else chrom
    ci = ref_names.index(chrom)
    ref_off = ref_offsets(ref_lengths)
    off = int(ref_off[ci])

    # sites annotes de la region (coordonnees locales 0-based)
    sp = np.asarray(ann.site_pos, dtype=np.int64) if ann is not None else np.zeros(0, np.int64)
    sel = (sp >= off + a1 - 1) & (sp <= off + b1)
    site_idx = np.where(sel)[0]
    sites_local = sp[site_idx] - off

    junc, nons, n = scan_region(bam, chrom, a1, b1, anchor=anchor,
                                min_mapq=min_mapq, sites_local=sites_local)
    log(f"{n} reads chevauchant {chrom}:{a1:,}-{b1:,}")

    # fichier de comptage de l'echantillon
    stored_j: Dict[int, int] = {}
    stored_ns = None
    cinfo = ""
    if sample:
        cp = wd.sample_counts(sample)
        if cp.exists():
            z = np.load(winpath(cp), allow_pickle=True)
            stored_j = dict(zip(z["jkey"].tolist(), z["jcount"].tolist()))
            stored_ns = z["nonsplit"]
            info = z["info"]
            cinfo = str(info[3]) if len(info) >= 4 else "(BAM non enregistre : comptage < 2.5)"
            log(f"fichier de comptage : {cp}")
            log(f"  BAM au moment du comptage : {cinfo}")
            from .count import bam_signature, _same_bam
            if len(info) >= 4 and str(info[3]) and not _same_bam(str(info[3]), bam_signature(bam)):
                warn("  LE BAM A CHANGE depuis le comptage : relancer count (le "
                     "comptage perime sera refait automatiquement en 2.5)")
        else:
            warn(f"pas de fichier de comptage pour {sample}")

    # matrice de la cohorte
    jindex = {}
    if st is not None:
        jindex = {int(k): i for i, k in enumerate(st.jkey.tolist())}
    si = st.samples.index(sample) if (st is not None and sample in st.samples) else None

    rows: List[dict] = []
    keys = sorted(set(junc.keys()) | {
        (int((k >> 24) - off), int(k & ((1 << 24) - 1))) for k in stored_j
        if off + a1 - 1 - 1_000_000 <= (k >> 24) <= off + b1
        and (k >> 24) + (k & ((1 << 24) - 1)) >= off + a1 - 1})
    for (s0, ln) in keys:
        if s0 + ln < a1 - 1 or s0 > b1:
            continue
        key = pack_junction(off + s0, ln)
        c = junc.get((s0, ln), {})
        col = jindex.get(int(key))
        mat = coh_med = coh_n = None
        if col is not None:
            Kc = np.asarray(st.Ks[:, col], dtype=np.int64)
            ref = st.reference.copy()
            if si is not None:
                mat = int(Kc[si])
                ref[si] = False
            coh_med = float(np.median(Kc[ref])) if ref.any() else None
            coh_n = int((Kc[ref] >= 3).sum())
        rows.append(dict(
            kind="junction", position=f"{chrom}:{s0 + 1}-{s0 + ln}",
            **{k: int(c.get(k, 0)) for k in CLASSES},
            stored=int(stored_j.get(int(key), 0)) if stored_j else None,
            matrix=mat, cohort_median=coh_med, cohort_n_ge3=coh_n,
            in_matrix=col is not None))
    for i, s in zip(site_idx.tolist(), sites_local.tolist()):
        c = nons.get(int(s), {})
        stv = int(stored_ns[i]) if stored_ns is not None and i < len(stored_ns) else None
        mat = coh_med = None
        if st is not None and st.Kn is not None and si is not None:
            Kn = np.asarray(st.Kn[:, i], dtype=np.int64)
            mat = int(Kn[si])
            ref = st.reference.copy()
            ref[si] = False
            coh_med = float(np.median(Kn[ref])) if ref.any() else None
        rows.append(dict(kind="unspliced_at_site", position=f"{chrom}:{int(s)}|{int(s) + 1}",
                         **{k: int(c.get(k, 0)) for k in CLASSES},
                         stored=stv, matrix=mat, cohort_median=coh_med,
                         cohort_n_ge3=None, in_matrix=True))
    # sites partages entre deux blocs (bug de regroupement < 2.5)
    if st is not None:
        cols = [jindex.get(int(pack_junction(off + int(r["position"].split(":")[1].split("-")[0]) - 1,
                                             int(r["position"].split("-")[1]) - int(r["position"].split(":")[1].split("-")[0]) + 1)))
                for r in rows if r["kind"] == "junction"]
        cols = [c for c in cols if c is not None]
        if cols:
            from .store import check_blocks
            c = np.array(cols)
            nb = check_blocks(st.jstart[c], st.jend[c], st.block_id[c])
            if nb:
                warn(f"{nb} site(s) de la region partage(s) entre deux blocs de la "
                     f"matrice : psi/theta FAUX ici (bug de regroupement des versions "
                     f"< 2.5). Relancer `build` puis `analyze` avec la 2.5.")
                for r in rows:
                    r["split_block_bug"] = True
    _print(rows)
    out = wd.results / f"probe_{safe_name(sample or Path(bam).stem)}_{chrom}_{a1}_{b1}.tsv"
    try:
        out.parent.mkdir(parents=True, exist_ok=True)
        cols = ["kind", "position"] + list(CLASSES) + ["stored", "matrix",
                                                       "cohort_median", "cohort_n_ge3", "in_matrix"]
        with open(winpath(out), "w", encoding="utf-8") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join("" if r.get(k) is None else str(r.get(k)) for k in cols) + "\n")
        log(f"tableau : {out}")
    except OSError:
        pass
    _verdict(rows)
    return rows


def _print(rows: List[dict]) -> None:
    hdr = (f"{'type':18s} {'position':28s} {'kept':>6s} {'dup':>5s} {'sec':>5s} "
           f"{'supp':>5s} {'mapq':>5s} {'stored':>7s} {'matrix':>7s} {'cohort_med':>10s}")
    log(hdr)
    for r in rows:
        if r["kind"] == "junction" and r["kept"] + r["duplicate"] + r["secondary"] < 2 \
                and not (r["stored"] or r["matrix"]):
            continue
        f = lambda v: "" if v is None else (f"{v:.0f}" if isinstance(v, float) else str(v))  # noqa: E731
        log(f"{'jonction' if r['kind'] == 'junction' else 'non coupes au site':18s} "
            f"{r['position']:28s} {r['kept']:6d} {r['duplicate']:5d} {r['secondary']:5d} "
            f"{r['supplementary']:5d} {r['low_mapq']:5d} {f(r['stored']):>7s} "
            f"{f(r['matrix']) if r['in_matrix'] else 'ABSENT':>7s} {f(r['cohort_median']):>10s}")


def _verdict(rows: List[dict]) -> None:
    bad = []
    if rows and rows[0].get("split_block_bug"):
        warn("VERDICT : matrice construite avec le bug de regroupement : relancer "
             "`build` puis `analyze` (pas besoin de recompter)")
        return
    for r in rows:
        if r["stored"] is not None and r["kept"] >= 5 and r["stored"] < 0.5 * r["kept"]:
            bad.append(f"{r['position']} : {r['kept']} reads dans le BAM, {r['stored']} "
                       f"dans le fichier de comptage")
        if r["kind"] == "junction" and r["kept"] >= 10 and not r["in_matrix"]:
            bad.append(f"{r['position']} : {r['kept']} reads mais jonction absente de "
                       f"la matrice (filtree a 'build' ?)")
        if r["kept"] + r["duplicate"] >= 10 and r["duplicate"] > r["kept"]:
            bad.append(f"{r['position']} : majorite de doublons PCR ({r['duplicate']} "
                       f"contre {r['kept']}) -> essayer --keep-duplicates")
    if bad:
        warn("VERDICT : le probleme est au COMPTAGE")
        for b in bad:
            warn("  " + b)
    else:
        log("VERDICT : BAM, comptage et matrice concordent ; si l'evenement manque, "
            "voir `why` (modele / seuils)")
