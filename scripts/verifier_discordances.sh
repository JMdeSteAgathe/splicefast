#!/usr/bin/env bash
# Verifie, dans un dossier de resultats SpliceFast, les evenements de reference
# (vus par FRASER2) qui manquaient ou etaient mal classes en 2.3.
#
#   scripts/verifier_discordances.sh /data/splicefast
#
# Pour chaque cas : la ligne la plus proche dans hits_all.tsv (coordonnees a
# +-2 pb), sa porte d'entree, son anomalie et son rang dans l'echantillon.
set -euo pipefail
W="${1:?usage: $0 DOSSIER_DE_TRAVAIL}"
python3 - "$W/results/hits_all.tsv" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
cases = [  # prefixe echantillon, gene, debut, fin (1-based, comme FRASER)
    ("ECH01", "TBCK",    106171271, 106193608),
    ("ECH02", "NIPBL",   37026328,  37027358),
    ("ECH03", "TUBGCP2", 133282343, 133283077),
    ("ECH04", "POLR2A",  7513244,   7513306),
    ("ECH05", "SMC1A",   53403894,  53405011),
    ("ECH06", "DYNC1H1", 102000145, 102000285),
]
def num(x, d=1.0):
    try: return float(x)
    except ValueError: return d
for pre, gene, s, e in cases:
    sub = [r for r in rows if r["sample"].startswith(pre)]
    if not sub:
        print(f"{pre} {gene}: echantillon absent des resultats"); continue
    order = sorted({r["anomaly_id"] for r in sub},
                   key=lambda a: min(num(r["anomaly_p"]) for r in sub if r["anomaly_id"] == a))
    hit = [r for r in sub if r["gene"] == gene and (
        (abs(int(r["start"]) - s) <= 2 and abs(int(r["end"] or 0) - e) <= 2)
        or abs(int(r["start"]) - s) <= 2 or abs(int(r["start"]) - (e + 1)) <= 2)]
    if not hit:
        print(f"{pre} {gene} {s}-{e}: ABSENT -> lancer : python3 -m splicefast why "
              f"-w DOSSIER --sample <nom exact> --gene {gene}")
        continue
    a = hit[0]["anomaly_id"]
    rank = order.index(a) + 1 if a in order else "?"
    print(f"{pre} {gene} {s}-{e}: TROUVE, anomalie {a} rang {rank}/{len(order)} "
          f"p_anomalie={hit[0]['anomaly_p']} preuves={hit[0]['evidence']}")
    for r in hit:
        print(f"    {r['metrics']:<22} {r['start']}-{r['end']} k={r['split_reads']} "
              f"n={r['total_reads']} delta={r['delta_psi']} p={r['p_value']} retenu_par={r['selected_by']}")
EOF
