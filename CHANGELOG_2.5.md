# SpliceFast 2.5

## 1. NIPBL (ECH02) : ce qu'est l'événement, et pourquoi il manquait

**L'événement.** L'extrait de BAM (viewer sashimi) le montre sans ambiguïté :

| | patient | S41 | S42 |
|---|---|---|---|
| jonction canonique exon 31 → 32 (chr5:37026328-37027358) | 63 | 140 | 160 |
| **jonction nouvelle chr5:37026485-37027358** | **32** | 0 | 0 |
| reads non coupés au donneur canonique (37026327) | 32 | 0 | 2 |

Le variant c.5808+159C>T (chr5:37026486) transforme la séquence intronique
`G C AAG` en `G T AAG` : un donneur consensus parfait à +157. L'exon 31 est
prolongé de 157 pb dans l'intron, ce qui correspond à la « rétention intronique
partielle » vue dans IGV. Cela fait trois signaux indépendants : la jonction nouvelle (ψ3 = 0,34
contre 0), la baisse du Jaccard de la jonction canonique (celle que FRASER a
rapportée), et la chute de theta au donneur canonique.

**Dans une cohorte simulée reproduisant exactement ce cas, la 2.4 comme la 2.5
le détectent** (ψ3, Jaccard : p = 1e-40 ; theta : 1e-20). Les trois signaux
manquaient pourtant ensemble sur votre run, ce qui désigne les données
stockées et non le modèle. Deux causes réelles ont été trouvées et corrigées :

1. **Regroupement inexact des jonctions (bug, toutes versions < 2.5).** Les
   jonctions qui partagent un site étaient regroupées par propagation
   d'étiquettes limitée à 200 itérations. Sur une longue chaîne de jonctions
   (grand gène, gènes voisins reliés par des jonctions de lecture continue,
   fréquents avec 355 échantillons), la propagation ne convergeait pas : à
   partir de 400 exons enchaînés, un même groupe était coupé en 3 ou 4 (test
   `tests/test_modele.py`). Si une frontière de bloc tombait entre les
   morceaux, ψ5/ψ3/Jaccard/theta de ces sites étaient calculés sur un groupe
   incomplet, et ce sans aucun message. Le regroupement est maintenant exact
   (composantes connexes), et `build` refuse de produire un store dont un
   site serait coupé. `analyze` signale un store construit par une ancienne
   version.
2. **Comptage périmé.** Le cache de comptage était identifié par le seul nom
   d'échantillon. Un BAM remplacé, ou un autre fichier portant le même nom,
   réutilisait silencieusement l'ancien comptage. Chaque comptage enregistre
   désormais le chemin, la taille et la date du BAM, et il est refait si le BAM a changé.

**Pour trancher sur votre run Linux** (une minute, sans relancer l'analyse) :

```bash
python3 -m splicefast probe -w $run_dir --sample ECH02-PAX_S61 \
    --region chr5:37025000-37028000
```

`probe` relit la région directement dans le BAM, ventile les reads (retenus,
doublons PCR, secondaires, MAPQ), les compare au fichier de comptage puis à la
matrice, et conclut : problème au comptage, store à reconstruire, ou rien à
signaler (auquel cas `why` prend le relais).

La 2.5 nomme aussi ce type d'événement : `exon_extension (partial intron
retention)`, avec le détail `exon extension: +157 bp into intron …`. Sur le
brin −, les libellés `novel_donor` et `novel_acceptor` étaient inversés : ils
sont maintenant corrects.

## 2. Événements récurrents (POLR2A…)

Pour chaque ligne, le rapport compte les échantillons de la cohorte qui portent le **même**
événement (p ≤ 1e-4, même sens, |Δψ| ≥ seuil/2). Les colonnes ajoutées sont
`n_carriers` et `other_carriers`. L'événement est marqué `recurrent = yes` à
partir de 4 porteurs (`--recurrent`). Il s'agit alors probablement d'un polymorphisme,
d'un artefact technique ou d'un effet du tissu. Le visualiseur les masque par défaut,
de même que les loci hypervariables HLA/IG/TR (`hypervariable_locus`).

## 3. Bruit par échantillon

`summary.tsv` contient un `noise_index` : le taux de tests très significatifs de
l'échantillon rapporté à la médiane de la cohorte (1 = typique, > 2 = élevé,
> 4 = très élevé), décliné pour les jonctions et pour la couverture. Le rapport cite les
échantillons très bruyants. Dans votre cohorte, les FIBRO, TISSU et certains
PAX à 1 400–1 500 événements devraient sans doute passer en `--test-only`.

## 4. PanelApp : tous les panels, avec la confiance

```bash
python3 -m splicefast panelapp --update -w $run_dir
```

Cette commande télécharge l'API PanelApp Australie (≈ 34 000 entrées gène × panel), puis
chaque ligne reçoit **tous** les panels de son gène, chacun avec sa confiance
(3 et 4 = vert, 2 = ambre, 1 = rouge) :
- `panels`, du type `Mendeliome [G]; Intellectual disability … [A]` ;
- `panels_green` ;
- `panelapp_best_confidence` ;
- `moi`, fusionné à partir des panels verts ;
- `phenotypes` ;
- `panel_details`, au format JSON, lu par le visualiseur.

Derrière un proxy hospitalier, utiliser `--insecure`, ou `scripts/panelapp_download.ps1` puis
`--from-json`. Les anciens exports `--panel` restent acceptés en secours.

## 5. Nouveaux résultats

- **`SpliceFast_viewer.html`** remplace `result.html`. On l'ouvre, on choisit le dossier
  `results`, puis on coche un ou plusieurs échantillons.
  - Rien n'est envoyé : les fichiers sont lus localement.
  - Une ligne par anomalie, de hauteur minimale. Le locus tient dans une seule colonne,
    avec un bouton copier.
  - Filtres : nom de panel (recherche), confiance, type, p, |Δψ|, reads, récurrence,
    HLA/IG/TR, référence. Ce dernier filtre fonctionne, contrairement à celui de la 2.4.
  - Un clic sur une ligne ouvre la structure du gène avec les arcs de jonction, la
    distribution de ψ dans la cohorte avec la valeur du patient et l'attendu, tous les
    signaux de l'anomalie, les porteurs, et les panels avec confiance, MOI et phénotypes.
  - Navigation au clavier ↑/↓, export TSV de la sélection.
  - 355 échantillons et 178 000 lignes se chargent en une dizaine de secondes.
- Les TSV sont en anglais. `<sample>.tsv` contient tout, y compris les colonnes JSON
  du visualiseur. `hits_all.tsv` contient la même chose sans le JSON. On trouve aussi
  `anomalies.tsv`, `summary.tsv`, `panels.tsv`, `run_info.json` et `results.xlsx`.

## 6. Windows

- Le paquet Windows contient un `sfscan.exe` compilé statiquement (htslib 1.21) :
  - comptage et couverture en une seule lecture des BAM, 10 à 30 fois plus rapide que
    le lecteur Python ;
  - résultats identiques bit à bit au lecteur Python (auto-test) ;
  - s'il ne se lance pas (antivirus), SpliceFast revient seul au lecteur Python.
- sfscan écrase désormais un fichier de sortie existant (sous Windows, `rename`
  échouait après un run interrompu).

## 7. Référence par motif de nom

`--reference-match PAX` (sur `build`, `all` ou `reference`) : seuls les
échantillons dont le nom contient le motif (expression régulière, sans tenir compte
de la casse) servent de norme. Le critère se combine avec `--min-reads 60e6`. Les autres
échantillons sont testés sans participer à la référence. Dans `all`, la référence est
fixée **avant** l'analyse, et l'analyse est refaite si la référence a changé.
