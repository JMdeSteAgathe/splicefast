"""
PanelApp (Australie) : tous les panels, tous les genes, avec leur niveau de
confiance.

Telechargement (une fois, puis a rafraichir de temps en temps) :

    python -m splicefast panelapp --update -w DOSSIER

interroge https://panelapp-aus.org/api/v1/genes/ page par page (~340 pages,
quelques minutes) et ecrit un cache compact `panelapp_genes.tsv.gz`. Le
rapport l'utilise ensuite automatiquement : chaque evenement recoit la liste
COMPLETE des panels ou figure son gene, avec la confiance de chacun.

Confiance PanelApp : 3 et 4 = vert, 2 = ambre, 1 = rouge, 0 = non classe.

Hors ligne / derriere un proxy capricieux : enregistrer les pages JSON de
l'API dans un dossier (page_1.json, page_2.json, ...) puis
    python -m splicefast panelapp --from-json DOSSIER_JSON -w DOSSIER
(script PowerShell fourni : scripts/panelapp_download.ps1).
"""
from __future__ import annotations

import gzip
import io
import json
import os
import re
import ssl
import time
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

from .util import log, warn, winpath

API = "https://panelapp-aus.org/api/v1/genes/"
CACHE_NAME = "panelapp_genes.tsv.gz"
FIELDS = ("symbol", "hgnc_id", "ensembl_id", "panel_id", "panel_name",
          "disease_group", "panel_version", "confidence", "moi", "phenotypes",
          "penetrance", "mode_of_pathogenicity", "tags")

CONF_LABEL = {"4": "green", "3": "green", "2": "amber", "1": "red", "0": "none"}
CONF_RANK = {"green": 0, "amber": 1, "red": 2, "none": 3}
CONF_SHORT = {"green": "G", "amber": "A", "red": "R", "none": "-"}


def conf_label(level) -> str:
    return CONF_LABEL.get(str(level).strip(), "none")


# --------------------------------------------------------------------------- #
# Telechargement
# --------------------------------------------------------------------------- #
def _get(url: str, timeout: float, insecure: bool, retries: int = 5) -> dict:
    ctx = ssl.create_default_context()
    if insecure:
        ctx.check_hostname = False
        ctx.verify_mode = ssl.CERT_NONE
    req = urllib.request.Request(url, headers={
        "accept": "application/json", "User-Agent": "splicefast-panelapp/2.5"})
    last = None
    for k in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=ctx) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            last = e
            if isinstance(e, ssl.SSLError) or "CERTIFICATE" in str(e).upper():
                raise SystemExit(
                    f"PanelApp : erreur de certificat SSL ({e}). Reseau "
                    f"hospitalier avec inspection SSL : relancer avec "
                    f"--insecure, ou utiliser scripts/panelapp_download.ps1 "
                    f"puis --from-json.")
            time.sleep(min(30, 2 ** k))
    raise SystemExit(f"PanelApp injoignable ({url}) : {last}. Verifier la "
                     f"connexion / le proxy, ou utiliser --from-json.")


def _rows_from_page(page: dict) -> Iterable[dict]:
    for g in page.get("results", []):
        gd = g.get("gene_data") or {}
        panel = g.get("panel") or {}
        ens = ""
        try:
            e38 = (gd.get("ensembl_genes") or {}).get("GRch38") or {}
            for _rel, v in e38.items():
                ens = v.get("ensembl_id", "") or ens
        except AttributeError:
            pass
        sym = gd.get("hgnc_symbol") or gd.get("gene_symbol") or g.get("entity_name") or ""
        if not sym:
            continue
        yield dict(
            symbol=sym, hgnc_id=gd.get("hgnc_id") or "", ensembl_id=ens,
            panel_id=str(panel.get("id", "")), panel_name=panel.get("name", ""),
            disease_group=panel.get("disease_group", "") or "",
            panel_version=str(panel.get("version", "")),
            confidence=str(g.get("confidence_level", "")),
            moi=g.get("mode_of_inheritance") or "",
            phenotypes=" | ".join(p.strip() for p in (g.get("phenotypes") or []) if p and p.strip()),
            penetrance=g.get("penetrance") or "",
            mode_of_pathogenicity=g.get("mode_of_pathogenicity") or "",
            tags=",".join(g.get("tags") or []))


def _clean(v) -> str:
    return str(v).replace("\t", " ").replace("\r", " ").replace("\n", " ")


def write_cache(rows: List[dict], out: Path, source: str) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    with gzip.open(winpath(tmp), "wt", encoding="utf-8") as fh:
        fh.write(f"# source={source} date={time.strftime('%Y-%m-%d')} n={len(rows)}\n")
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(_clean(r.get(f, "")) for f in FIELDS) + "\n")
    os.replace(winpath(tmp), winpath(out))
    np_ = len({r["panel_id"] for r in rows})
    ng = len({r["symbol"] for r in rows})
    log(f"PanelApp : {len(rows)} entrees, {ng} genes, {np_} panels -> {out}")


def download(out: Path, api: str = API, timeout: float = 60.0,
             insecure: bool = False, save_json: Optional[Path] = None) -> Path:
    url = api
    rows: List[dict] = []
    n = 0
    total = None
    t0 = time.time()
    while url:
        page = _get(url, timeout, insecure)
        n += 1
        if total is None:
            total = page.get("count")
            log(f"PanelApp : {total} entrees a telecharger")
        if save_json is not None:
            save_json.mkdir(parents=True, exist_ok=True)
            (save_json / f"page_{n}.json").write_text(json.dumps(page), encoding="utf-8")
        rows.extend(_rows_from_page(page))
        url = page.get("next")
        if n % 20 == 0:
            log(f"  page {n} : {len(rows)} entrees ({time.time() - t0:.0f} s)")
    write_cache(rows, out, api)
    return out


def from_json(folder: Path, out: Path) -> Path:
    files = sorted(Path(folder).glob("*.json"),
                   key=lambda p: [int(x) if x.isdigit() else x
                                  for x in re.split(r"(\d+)", p.name)])
    if not files:
        raise SystemExit(f"aucun fichier .json dans {folder}")
    rows: List[dict] = []
    for f in files:
        rows.extend(_rows_from_page(json.loads(f.read_text(encoding="utf-8-sig"))))
    write_cache(rows, out, f"json:{folder}")
    return out


# --------------------------------------------------------------------------- #
# Annotation
# --------------------------------------------------------------------------- #
_MOI_ATOMS = [
    (r"^BOTH monoallelic and biallelic", ("AD", "AR")),
    (r"^BIALLELIC", ("AR",)),
    (r"^MONOALLELIC", ("AD",)),
    (r"^X[- ]?LINKED.*monoallelic mutations in females", ("XLD",)),
    (r"^X[- ]?LINKED", ("XLR",)),
    (r"^MITOCHONDRIAL", ("MT",)),
    (r"^Other", ("other",)),
    (r"^Unknown", ("?",)),
]
_MOI_ORDER = ["AD", "AR", "XLD", "XLR", "MT", "other", "?"]


def moi_short(text: str) -> str:
    t = (text or "").strip()
    if not t:
        return ""
    for pat, atoms in _MOI_ATOMS:
        if re.match(pat, t, flags=re.I):
            s = "/".join(atoms)
            if "imprint" in t.lower() and "unknown" not in t.lower():
                s += " (imprinted)"
            return s
    return "other"


def merge_moi(texts: Iterable[str]) -> str:
    got = set()
    for t in texts:
        for pat, atoms in _MOI_ATOMS:
            if re.match(pat, (t or "").strip(), flags=re.I):
                got.update(atoms)
                break
    return "/".join([a for a in _MOI_ORDER if a in got])


class PanelAppDB:
    """Gene -> panels (avec confiance), depuis le cache telecharge."""

    def __init__(self) -> None:
        self.by_symbol: Dict[str, List[dict]] = defaultdict(list)
        self.by_ensembl: Dict[str, List[dict]] = defaultdict(list)
        self.panels: Dict[str, dict] = {}
        self.source = ""

    @classmethod
    def load(cls, path: Path) -> "PanelAppDB":
        db = cls()
        with gzip.open(winpath(path), "rt", encoding="utf-8") as fh:
            first = fh.readline()
            if first.startswith("#"):
                db.source = first[1:].strip()
                header = fh.readline().rstrip("\n").split("\t")
            else:
                header = first.rstrip("\n").split("\t")
            for line in fh:
                v = line.rstrip("\n").split("\t")
                r = dict(zip(header, v))
                r["conf"] = conf_label(r.get("confidence", ""))
                db.by_symbol[r["symbol"].upper()].append(r)
                if r.get("ensembl_id"):
                    db.by_ensembl[r["ensembl_id"].split(".")[0]].append(r)
                p = db.panels.setdefault(r["panel_id"], dict(
                    panel_id=r["panel_id"], panel_name=r["panel_name"],
                    disease_group=r.get("disease_group", ""),
                    version=r.get("panel_version", ""), n_green=0, n_genes=0))
                p["n_genes"] += 1
                p["n_green"] += r["conf"] == "green"
        log(f"PanelApp : {sum(len(v) for v in db.by_symbol.values())} entrees, "
            f"{len(db.by_symbol)} genes, {len(db.panels)} panels ({db.source})")
        return db

    def entries(self, symbol: str, gene_id: str = "") -> List[dict]:
        e = self.by_symbol.get((symbol or "").upper())
        if not e and gene_id:
            e = self.by_ensembl.get(gene_id.split(".")[0])
        return e or []

    def annotate(self, symbol: str, gene_id: str = "") -> dict:
        ent = sorted(self.entries(symbol, gene_id),
                     key=lambda r: (CONF_RANK[r["conf"]], r["panel_name"].lower()))
        if not ent:
            return dict(panels="", panels_green="", panel_confidence="",
                        n_panels=0, moi="", phenotypes="", panel_details="")
        best = ent[0]["conf"]
        green = [r for r in ent if r["conf"] == "green"]
        use = green or [r for r in ent if r["conf"] == "amber"] or ent
        phen: List[str] = []
        for r in use:
            for ph in r.get("phenotypes", "").split(" | "):
                ph = ph.strip()
                if ph and ph not in phen:
                    phen.append(ph)
        details = [[r["panel_name"], CONF_SHORT[r["conf"]], moi_short(r.get("moi", "")),
                    r.get("phenotypes", "")] for r in ent]
        return dict(
            panels="; ".join(f"{r['panel_name']} [{CONF_SHORT[r['conf']]}]" for r in ent),
            panels_green="; ".join(r["panel_name"] for r in green),
            panel_confidence=best if best != "none" else "",
            n_panels=len(ent),
            moi=merge_moi(r.get("moi", "") for r in use),
            phenotypes=" | ".join(phen),
            panel_details=json.dumps(details, ensure_ascii=False, separators=(",", ":")))

    def write_panel_list(self, out: Path) -> None:
        with open(winpath(out), "w", encoding="utf-8") as fh:
            fh.write("panel_id\tpanel_name\tdisease_group\tversion\tn_genes\tn_green\n")
            for p in sorted(self.panels.values(), key=lambda x: x["panel_name"].lower()):
                fh.write("\t".join(_clean(p[k]) for k in ("panel_id", "panel_name", "disease_group",
                                                          "version", "n_genes", "n_green")) + "\n")


def default_cache(wd) -> Path:
    return Path(wd.root) / "panelapp" / CACHE_NAME


def find_cache(wd, explicit: Optional[str] = None) -> Optional[Path]:
    """Cache explicite, sinon celui du dossier de travail, sinon ~/.splicefast."""
    cands = [Path(explicit)] if explicit else []
    cands += [default_cache(wd), Path.home() / ".splicefast" / CACHE_NAME]
    for c in cands:
        if c.exists():
            return c
    if explicit:
        warn(f"cache PanelApp introuvable : {explicit}")
    return None
