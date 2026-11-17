#!/usr/bin/env bash
# Verifie, dans un dossier de resultats SpliceFast, les evenements de reference
# (vus par FRASER2) qui manquaient ou etaient mal classes en 2.3.
#
#   scripts/verifier_discordances.sh /data/splicefast [CAS.tsv]
#
# CAS.tsv : une ligne par cas, 4 colonnes separees par des tabulations :
#   prefixe_echantillon  gene  debut  fin      (1-based, comme FRASER)
# Par defaut : scripts/discordances.local.tsv. Ce fichier contient de vrais
# identifiants d'echantillons : il n'est PAS versionne (.gitignore). Modele :
# scripts/discordances.exemple.tsv.
#
# Pour chaque cas : la ligne la plus proche dans hits_all.tsv (coordonnees a
# +-2 pb), sa porte d'entree, son anomalie et son rang dans l'echantillon.
set -euo pipefail
W="${1:?usage: $0 DOSSIER_DE_TRAVAIL [CAS.tsv]}"
CAS="${2:-$(dirname "$0")/discordances.local.tsv}"
[ -f "$CAS" ] || { echo "fichier de cas introuvable : $CAS (modele : scripts/discordances.exemple.tsv)" >&2; exit 1; }
python3 - "$W/results/hits_all.tsv" "$CAS" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1], encoding="utf-8"), delimiter="\t"))
cases = []  # prefixe echantillon, gene, debut, fin (1-based, comme FRASER)
for line in open(sys.argv[2], encoding="utf-8"):
    f = line.rstrip("\r\n").split("\t")
    if not line.strip() or line.startswith("#") or len(f) < 4 or not f[2].strip().isdigit():
        continue  # commentaire, ligne vide ou en-tete
    cases.append((f[0].strip(), f[1].strip(), int(f[2]), int(f[3])))
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
