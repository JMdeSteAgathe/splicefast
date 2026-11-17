"""
Annotation des genes par les panels PanelApp (Incidentalome, Mendeliome, ...).

Lit les exports TSV de PanelApp et fournit, pour chaque symbole de gene :
    * la liste des panels ou il figure,
    * le(s) mode(s) de transmission (Model_Of_Inheritance),
    * les phenotypes associes (fusionnes, dedoublonnes),
    * le niveau de confiance GEL (3 = vert, 2 = ambre, 1 = rouge).

Un gene present dans plusieurs panels n'apparait qu'une fois : les panels sont
cumules et les phenotypes fusionnes. Si les modes de transmission different
entre panels, ils sont tous conserves et prefixes par le nom du panel.
"""
from __future__ import annotations

import csv
import io
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .util import log, warn, winpath

# Abreviations lisibles des libelles PanelApp (qui sont tres verbeux).
_MOI_SHORT = [
    (r"^BIALLELIC", "AR"),
    (r"^BOTH monoallelic and biallelic.*imprinted", "AD/AR"),
    (r"^BOTH monoallelic and biallelic", "AD/AR"),
    (r"^MONOALLELIC.*imprinted, paternally", "AD (imprint. pat.)"),
    (r"^MONOALLELIC.*imprinted, maternally", "AD (imprint. mat.)"),
    (r"^MONOALLELIC", "AD"),
    (r"^X[- ]?LINKED.*monoallelic mutations in females", "XLD"),
    (r"^X[- ]?LINKED", "XLR"),
    (r"^MITOCHONDRIAL", "MITO"),
    (r"^Other.*", "autre"),
    (r"^Unknown", "?"),
]

_GEL_LABEL = {"3": "vert", "2": "ambre", "1": "rouge", "0": "?"}

# Decomposition d'un libelle PanelApp en modes elementaires, pour pouvoir
# fusionner proprement plusieurs panels : "BOTH mono and biallelic" vaut
# {AD, AR}, donc "both" + "AD" -> AD/AR, et "AD" + "AD" -> AD.
_MOI_ATOMS = [
    (r"^BOTH monoallelic and biallelic", ("AD", "AR")),
    (r"^BIALLELIC", ("AR",)),
    (r"^MONOALLELIC", ("AD",)),
    (r"^X[- ]?LINKED.*monoallelic mutations in females", ("XLD",)),
    (r"^X[- ]?LINKED", ("XLR",)),
    (r"^MITOCHONDRIAL", ("MITO",)),
    (r"^Other", ("autre",)),
    (r"^Unknown", ("?",)),
]
_MOI_ORDER = ["AD", "AR", "XLD", "XLR", "MITO", "autre", "?"]


def moi_atoms(text: str) -> tuple:
    """Modes elementaires d'un libelle PanelApp (ex: -> ('AD', 'AR'))."""
    t = (text or "").strip()
    if not t:
        return ()
    for pat, atoms in _MOI_ATOMS:
        if re.match(pat, t, flags=re.I):
            return atoms
    return ("autre",)


def merge_moi(texts) -> str:
    """Fusionne plusieurs libelles en une transmission unique et lisible."""
    got = set()
    for t in texts:
        got.update(moi_atoms(t))
    if not got:
        return ""
    out = [a for a in _MOI_ORDER if a in got]
    out += sorted(got - set(_MOI_ORDER))
    return "/".join(out)


def short_moi(text: str) -> str:
    """Abrege un libelle de transmission PanelApp (ex: -> 'AD', 'AR', 'XLR')."""
    t = (text or "").strip()
    if not t:
        return ""
    for pat, short in _MOI_SHORT:
        if re.match(pat, t, flags=re.I):
            return short
    return t[:24]


def _split_list(text: str) -> List[str]:
    out = []
    for part in re.split(r"[;|]", text or ""):
        p = part.strip()
        if p and p not in out:
            out.append(p)
    return out


class PanelSet:
    """Annotation genique issue d'un ou plusieurs panels PanelApp."""

    def __init__(self) -> None:
        self.panels: List[str] = []
        self.genes: Dict[str, dict] = {}

    # ------------------------------------------------------------------ #
    def load(self, path: str, name: Optional[str] = None) -> int:
        p = Path(path)
        panel = name or p.stem
        try:
            raw = open(winpath(str(p)), "rb").read()
        except OSError as e:
            warn(f"panel illisible ({path}) : {e}")
            return 0
        text = raw.decode("utf-8-sig", errors="replace").replace("\r\n", "\n")
        rd = csv.DictReader(io.StringIO(text), delimiter="\t")
        if not rd.fieldnames:
            warn(f"panel vide : {path}")
            return 0
        cols = {c.strip().lower(): c for c in rd.fieldnames}

        def col(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        c_sym = col("gene symbol", "entity name", "gene_symbol")
        c_moi = col("model_of_inheritance", "model of inheritance",
                    "mode of inheritance")
        c_phe = col("phenotypes", "phenotype")
        c_gel = col("gel_status", "gel status")
        c_omim = col("omim")
        c_hpo = col("hpo")
        if c_sym is None:
            warn(f"colonne 'Gene Symbol' absente de {path}")
            return 0
        n = 0
        for row in rd:
            sym = (row.get(c_sym) or "").strip().upper()
            if not sym:
                continue
            g = self.genes.setdefault(sym, dict(panels=[], moi={}, phen={},
                                                gel={}, omim=[], hpo=[]))
            if panel not in g["panels"]:
                g["panels"].append(panel)
            moi = (row.get(c_moi) or "").strip() if c_moi else ""
            if moi:
                g["moi"].setdefault(panel, moi)
            if c_phe:
                lst = g["phen"].setdefault(panel, [])
                for ph in _split_list(row.get(c_phe) or ""):
                    if ph not in lst:
                        lst.append(ph)
            if c_gel:
                gel = (row.get(c_gel) or "").strip()
                if gel:
                    g["gel"][panel] = gel
            if c_omim:
                for o in _split_list(row.get(c_omim) or ""):
                    if o not in g["omim"]:
                        g["omim"].append(o)
            if c_hpo:
                for o in _split_list(row.get(c_hpo) or ""):
                    if o not in g["hpo"]:
                        g["hpo"].append(o)
            n += 1
        if panel not in self.panels:
            self.panels.append(panel)
        log(f"panel {panel} : {n} entrees, {len(self.genes)} genes au total")
        return n

    # ------------------------------------------------------------------ #
    def _abbrev(self) -> Dict[str, str]:
        """Abreviations courtes et NON ambigues des noms de panels."""
        if getattr(self, "_abbr", None) and self._abbr_for == tuple(self.panels):
            return self._abbr
        for k in range(3, 20):
            cand = {p: p[:k] for p in self.panels}
            if len(set(cand.values())) == len(cand):
                break
        self._abbr = cand
        self._abbr_for = tuple(self.panels)
        return cand

    def annotate(self, gene: str) -> dict:
        """Renvoie les champs prets a afficher pour un symbole de gene."""
        g = self.genes.get((gene or "").strip().upper())
        if not g:
            return dict(panels="", moi="", moi_long="", phenotypes="",
                        panel_conf="", omim="", hpo="", source_panel="")
        # PRIORITE : le dernier panel passe en ligne de commande gagne.
        # Pour un gene present dans plusieurs panels, transmission et
        # phenotypes viennent du panel le plus prioritaire qui les renseigne
        # (pas de melange : les libelles des deux panels ne concordent pas
        # toujours, et melanger produit des lignes illisibles).
        order = [p for p in reversed(self.panels) if p in g["panels"]]
        src_moi = next((p for p in order if g["moi"].get(p)), None)
        moi_s = merge_moi([g["moi"][src_moi]]) if src_moi else ""
        moi_l = g["moi"].get(src_moi, "") if src_moi else ""
        src_phe = next((p for p in order if g["phen"].get(p)), None)
        phen = g["phen"].get(src_phe, []) if src_phe else []
        conf = list(dict.fromkeys(_GEL_LABEL.get(v, v) for v in g["gel"].values()))
        return dict(panels=",".join(g["panels"]),
                    moi=moi_s, moi_long=moi_l,
                    phenotypes="; ".join(phen),
                    source_panel=src_moi or src_phe or "",
                    panel_conf="/".join(conf),
                    omim=";".join(g["omim"]), hpo=";".join(g["hpo"]))

    def __len__(self) -> int:
        return len(self.genes)


def load_panels(paths: Sequence[str]) -> Optional[PanelSet]:
    """Charge plusieurs fichiers de panel. Renvoie None si aucun n'est lisible."""
    if not paths:
        return None
    ps = PanelSet()
    for p in paths:
        if not p:
            continue
        if not os.path.exists(winpath(p)):
            warn(f"panel introuvable, ignore : {p}")
            continue
        ps.load(p)
    return ps if len(ps) else None


# --------------------------------------------------------------------------- #
# Identification du run d'un echantillon a partir du chemin de son BAM
# --------------------------------------------------------------------------- #
_RUN_RE = re.compile(r"(?i)^.*?((?:novaseq|nextseq|miseq|hiseq|run)[\w.\-]*)$")
# dossiers trop generiques pour servir de nom de run
_GENERIC = {"bam", "bams", "align", "aligned", "alignment", "alignments",
            "star", "output", "outputs", "results", "data", "mapping", "cram"}


def infer_run(bam_path: str) -> str:
    """Devine le nom du run depuis le chemin du BAM.

    Prend le dernier dossier du chemin qui ressemble a un run
    (ex: .../RNA-seq_total/novaseq_run4/xxx/sample.bam -> 'novaseq_run4'),
    sinon le dossier parent du BAM.
    """
    if not bam_path:
        return ""
    parts = [q for q in re.split(r"[\\/]+", str(bam_path)) if q]
    if not parts:
        return ""
    for q in reversed(parts[:-1]):
        m = _RUN_RE.match(q)
        if m:
            return m.group(1)
    for q in reversed(parts[:-1]):
        if q.lower().strip("._-") not in _GENERIC and not q.endswith(":"):
            return q
    return ""
