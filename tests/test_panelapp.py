"""Lecture d'une page de l'API PanelApp Australie et annotation d'un gene.

    python3 tests/test_panelapp.py
"""
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from splicefast.panelapp import PanelAppDB, from_json  # noqa: E402


def entry(sym, pid, name, conf, moi, phen):
    return {"gene_data": {"hgnc_symbol": sym, "hgnc_id": "HGNC:1",
                          "ensembl_genes": {"GRch38": {"90": {"ensembl_id": "ENSG00000164190"}}}},
            "entity_type": "gene", "entity_name": sym, "confidence_level": conf,
            "phenotypes": phen, "mode_of_inheritance": moi, "tags": [],
            "panel": {"id": pid, "name": name, "disease_group": "", "version": "1.0"}}


page = {"count": 3, "next": None, "results": [
    entry("NIPBL", 1, "Cornelia de Lange syndrome", "3",
          "MONOALLELIC, autosomal or pseudoautosomal", ["Cornelia de Lange syndrome 1, MIM# 122470"]),
    entry("NIPBL", 2, "Mendeliome", "4", "MONOALLELIC, autosomal or pseudoautosomal", ["CdLS"]),
    entry("NIPBL", 3, "Some screening panel", "2", "BIALLELIC, autosomal or pseudoautosomal", []),
]}
with tempfile.TemporaryDirectory() as td:
    (Path(td) / "page_1.json").write_text(json.dumps(page))
    out = Path(td) / "c.tsv.gz"
    from_json(Path(td), out)
    db = PanelAppDB.load(out)
    a = db.annotate("NIPBL")
    b = db.annotate("XXX", "ENSG00000164190.20")   # par identifiant Ensembl
    ok = (a["panel_confidence"] == "green" and a["panels_green"].count(";") == 1
          and "[A]" in a["panels"] and a["moi"] == "AD" and b["n_panels"] == 3)
    print(a["panels"], "|", a["moi"], "|", a["phenotypes"])
    print("test panelapp :", "OK" if ok else "ECHEC")
    sys.exit(0 if ok else 1)
