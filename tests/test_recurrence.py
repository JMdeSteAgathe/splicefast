"""2.7 : drapeau hypervariable par nom de gene, et porteurs plus stricts (cas SYNGAP1)."""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from splicefast.report import _carrier_idx, _set_carriers, hypervariable  # noqa: E402

# 1. Region CMH (chr6:28-34 Mb) : seuls HLA/IG/TR/KIR sont hypervariables
for g in ["SYNGAP1", "DAXX", "TAP1", "COL11A2", "TRAF3", "KIRREL1", "IGF1", "TRAPPC9"]:
    assert not hypervariable({"gene": g}), g
for g in ["HLA-A", "HLA-DRB1", "IGKC", "IGHG1", "IGHV3-23", "TRAC", "TRBV20-1", "KIR2DL1"]:
    assert hypervariable({"gene": g}), g
print("[1] hypervariable : SYNGAP1 non marque, HLA/IG/TR/KIR marques -> OK")

# 2. Porteurs : un echantillon a p 1e-3 ou de Delta faible n'est plus compte ;
#    un echantillon bruyant est liste a part, sans compter.
names = np.array(["P", "A", "B", "C", "N"])
mlog = np.array([30.0, 3.0, 12.0, 8.0, 20.0])        # p patient 1e-30, A 1e-3, ...
delta = np.array([-0.75, -0.70, -0.20, -0.50, -0.60])  # B : Delta trop faible (< 0.375)
ok = np.ones(5, bool)
h = {"sample": "P", "delta": -0.75}
car = _carrier_idx(h, mlog, delta, ok, lim=4.0, delta_cutoff=0.1, self_idx=0)
assert [names[i] for i in car] == ["P", "C", "N"], car
_set_carriers(h, car, names, noisy={"N"}, ref=np.array([1, 1, 1, 1, 1], bool),
              getkn=lambda i: (10, 100), mlog10p=mlog, delta=delta)
assert h["n_carriers"] == 2 and h["n_carriers_noisy"] == 1, h
assert "C" in h["carriers"] and "[+1 noisy]" in h["carriers"], h["carriers"]
print(f"[2] porteurs : {h['n_carriers']} (+{h['n_carriers_noisy']} bruyant) '{h['carriers']}' -> OK")
print("test recurrence : OK")
