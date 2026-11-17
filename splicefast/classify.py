"""
Etape 4a : classification biologique des evenements significatifs.

Types produits :
  exon_skipping        jonction dont les 2 sites sont annotes mais qui saute
                       au moins un exon annote (intron non annote)
  pseudo_exon          paire de jonctions nouvelles delimitant un exon
                       cryptique a l'interieur d'un intron annote
  novel_donor          nouveau site donneur (5')
  novel_acceptor       nouveau site accepteur (3')
  splice_site_shift    nouveau site a moins de N pb d'un site annote
                       (decalage de site, souvent variant +/-1,2)
  novel_junction       les deux sites sont nouveaux
  intron_retention     chute de theta = plus de reads a cheval exon/intron
  efficiency_gain      hausse de theta (moins de retention que la cohorte)
  annotated_intron     intron annote dont l'usage change (exon skipping deja
                       annote, isoforme alternative...)
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from .gtf import Annotation


class Classifier:
    def __init__(self, ann: Annotation, ref_names, ref_off,
                 shift_window: int = 20, pseudo_exon_max: int = 1000,
                 pseudo_exon_min: int = 15):
        self.ann = ann
        self.ref_names = list(ref_names)
        self.ref_off = np.asarray(ref_off)
        self.shift_window = shift_window
        self.pex_max = pseudo_exon_max
        self.pex_min = pseudo_exon_min
        self.has_ann = ann.n_sites > 0
        if self.has_ann:
            o = np.argsort(ann.exon_start)
            self.ex_s = ann.exon_start[o]
            self.ex_e = ann.exon_end[o]
            self.ex_g = ann.exon_gene[o]
            self.ex_emax = np.maximum.accumulate(self.ex_e)
            io = np.argsort(ann.intron_start)
            self.in_s = ann.intron_start[io]
            self.in_e = ann.intron_end[io]
            self.in_g = ann.intron_gene[io]
            self.in_emax = np.maximum.accumulate(self.in_e)

    # ------------------------------------------------------------------ #
    def gene_label(self, gi: int) -> Tuple[str, str, str]:
        if gi is None or gi < 0 or gi >= len(self.ann.gene_name):
            return ("", "", ".")
        strand = self.ann.gene_strand[gi]
        return (self.ann.gene_name[gi], self.ann.gene_id[gi],
                "+" if strand > 0 else ("-" if strand < 0 else "."))

    def local(self, gpos: int) -> Tuple[str, int]:
        ci = int(np.searchsorted(self.ref_off, gpos, side="right") - 1)
        return self.ref_names[ci], int(gpos - self.ref_off[ci])

    def dist_to_site(self, gpos: int) -> int:
        if not self.has_ann:
            return -1
        sp = self.ann.site_pos
        i = int(np.searchsorted(sp, gpos))
        best = 10 ** 9
        for j in (i - 1, i, i + 1):
            if 0 <= j < len(sp):
                best = min(best, abs(int(sp[j]) - gpos))
        return best

    def n_skipped_exons(self, gstart: int, gend: int, gene: int) -> int:
        if not self.has_ann:
            return 0
        i0 = int(np.searchsorted(self.ex_s, gstart, side="left"))
        i1 = int(np.searchsorted(self.ex_s, gend, side="right"))
        if i1 <= i0:
            return 0
        sel = slice(i0, i1)
        inside = (self.ex_e[sel] <= gend) & (self.ex_s[sel] >= gstart)
        if gene is not None and gene >= 0:
            inside &= (self.ex_g[sel] == gene)
        return int(inside.sum())

    def enclosing_intron(self, gpos: int) -> int:
        """Indice de l'intron annote contenant strictement gpos, sinon -1."""
        if not self.has_ann:
            return -1
        i = int(np.searchsorted(self.in_s, gpos, side="right") - 1)
        steps = 0
        while i >= 0 and self.in_emax[i] > gpos and steps < 500:
            if self.in_s[i] < gpos < self.in_e[i]:
                return i
            i -= 1
            steps += 1
        return -1

    # ------------------------------------------------------------------ #
    def classify_junction(self, gstart: int, gend: int, don_known: bool,
                          acc_known: bool, known_intron: bool,
                          gene: int) -> Dict:
        out: Dict[str, object] = {}
        nsk = 0
        if known_intron:
            etype = "annotated_intron"
            nsk = self.n_skipped_exons(gstart, gend, gene)
        elif don_known and acc_known:
            nsk = self.n_skipped_exons(gstart, gend, gene)
            etype = "exon_skipping" if nsk >= 1 else "novel_combination"
        elif don_known and not acc_known:
            d = self.dist_to_site(gend)
            etype = "splice_site_shift" if 0 < d <= self.shift_window else "novel_acceptor"
            out["dist_annot"] = d
        elif acc_known and not don_known:
            d = self.dist_to_site(gstart)
            etype = "splice_site_shift" if 0 < d <= self.shift_window else "novel_donor"
            out["dist_annot"] = d
        else:
            etype = "novel_junction"
        out["type"] = etype
        out["n_skipped_exons"] = nsk
        return out

    def pair_pseudo_exons(self, hits: List[dict]) -> None:
        """Detecte les paires de jonctions nouvelles encadrant un pseudo-exon.

        Cas typique : donneur annote -> accepteur nouveau A, puis donneur
        nouveau D -> accepteur annote, avec A < D dans le meme intron.
        L'exon cryptique est alors [A, D).
        """
        left: List[Tuple[int, dict]] = []   # (A, hit) accepteur nouveau
        right: List[Tuple[int, dict]] = []  # (D, hit) donneur nouveau
        for h in hits:
            if h.get("metric") == "theta" or h.get("gend") is None:
                continue
            if not h.get("acc_known", True) and h.get("don_known", False):
                left.append((int(h["gend"]), h))
            if not h.get("don_known", True) and h.get("acc_known", False):
                right.append((int(h["gstart"]), h))
            if not h.get("don_known", True) and not h.get("acc_known", True):
                left.append((int(h["gend"]), h))
                right.append((int(h["gstart"]), h))
        for a, ha in left:
            for d, hd in right:
                if ha is hd:
                    continue
                L = d - a
                if not (self.pex_min <= L <= self.pex_max):
                    continue
                i1 = self.enclosing_intron(a)
                i2 = self.enclosing_intron(d)
                if i1 < 0 or i1 != i2:
                    continue
                chrom, la = self.local(a)
                _, ld = self.local(d)
                note = f"pseudo-exon {chrom}:{la+1}-{ld} ({L} pb)"
                for h in (ha, hd):
                    h["type"] = "pseudo_exon"
                    h["comment"] = (h.get("comment", "") + "; " + note).strip("; ")
