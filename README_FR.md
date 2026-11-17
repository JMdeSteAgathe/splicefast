# SpliceFast 2.7 — Linux et Windows, lecteur C++

Détection d'événements d'épissage aberrants (saut d'exon, pseudo-exon, site
cryptique, extension d'exon, **rétention d'intron**) dans un run RNA-seq, par
comparaison de chaque individu à toute la cohorte — la même analyse que FRASER,
sans jamais saturer la mémoire.

**2.7** (`CHANGELOG_2.7.md`) : les gènes de la région du CMH (SYNGAP1…) ne sont
plus marqués « hypervariables », et la récurrence est plus stricte (les
échantillons bruyants ne comptent plus). Le visualiseur affiche la qualité de
l'échantillon et un tableau des porteurs. Relancer `report` seulement.

**2.6 — mise à jour indispensable** (`CHANGELOG_2.6.md`) : une partie des
comptages n'était pas placée dans la matrice (jonctions canoniques à 0 read, d'où
de faux ψ = 1). Relancer `build`, `analyze` et `report`, sans recompter.

**Nouveau en 2.5** (détail : `CHANGELOG_2.5.md`) :
- **visualiseur** `results/SpliceFast_viewer.html` : il charge le dossier de résultats et
  permet d'interroger un ou plusieurs échantillons, avec le détail et des mini-graphiques
  au clic ;
- **PanelApp** : tous les panels de l'API australienne, avec la confiance de chaque gène
  (`splicefast panelapp --update`) ;
- **récurrence** : nombre d'échantillons porteurs du même événement, le visualiseur masque
  les événements récurrents par défaut ;
- **indice de bruit** par échantillon dans `summary.tsv` ;
- **`probe`** : relit une région dans le BAM et la compare au comptage stocké ;
- **regroupement des jonctions exact** (bug des versions < 2.5) : relancer `build`
  puis `analyze`, sans recompter ;
- **Windows** : `sfscan.exe` fourni (installation : `INSTALL_WINDOWS.md`).

---

## 1. Pourquoi FRASER plante et pourquoi SpliceFast ne plante pas

FRASER charge en mémoire une dizaine de matrices `double` de taille
*jonctions × échantillons* (K, N pour ψ5, ψ3, θ, les valeurs attendues, les
poids de l'auto-encodeur…). Avec 3 millions de jonctions et 80 échantillons,
chaque matrice fait déjà **2 Go**, et l'auto-encodeur en manipule plusieurs
simultanément : on dépasse les 30–50 Go.

SpliceFast prend le problème dans l'autre sens :

| | FRASER | SpliceFast |
|---|---|---|
| Stockage | HDF5, `double`, jonctions × échantillons | memmap `int32`, **échantillons × jonctions** |
| Unité de calcul | toute la matrice | un **bloc** de ~20 000 jonctions |
| RAM pendant l'analyse | 30–50 Go | **~150 Mo**, indépendamment du nombre de jonctions |
| Modèle | auto-encodeur (itératif, toute la matrice) | SVD + régression robuste par jonction |
| Dépendances | R + Bioconductor + Rsubread | Python + numpy + scipy (+ htslib pour le lecteur C++) |

Le découpage en blocs est **exact et non approximatif** : les jonctions sont
d'abord regroupées en *clusters* (composantes connexes partageant un site
donneur ou accepteur), et un cluster n'est jamais coupé entre deux blocs. Les
dénominateurs de ψ5/ψ3/Jaccard sont donc toujours calculés sur la totalité du
groupe concerné, exactement comme si tout était en mémoire.

C'est votre intuition de départ (« gène par gène, reformater la matrice en
amont ») poussée à sa forme la plus efficace : le regroupement se fait par
cluster d'épissage plutôt que par gène (plus petit, plus sûr, et il gère les
jonctions inter-gènes), et la matrice est transposée à l'écriture pour que
lire un bloc coûte S petites lectures contiguës.

---

## 2. Installation (Linux, sans droits root)

Il faut Python 3.8+, numpy et scipy, et un compilateur C++ (g++ ≥ 7) :

```bash
python3 -m pip install --user numpy scipy
```

Décompressez l'archive, par exemple dans `~/splicefast`, puis compilez le
lecteur C++. Trois possibilités, de la plus simple à la plus autonome :

```bash
cd ~/splicefast

# a) htslib installée par le système ou conda
sudo apt install libhts-dev        # ou : conda install -c bioconda htslib
make

# b) aucun droit : compile htslib (+ libdeflate si cmake est présent) dans ./third_party
./install_htslib.sh

# c) htslib déjà compilée ailleurs
make HTSLIB=/chemin/vers/htslib-1.21
```

`libdeflate` divise environ par deux le temps de décompression ; si `cmake`
manque, `python3 -m pip install --user cmake` avant `./install_htslib.sh`.
Le binaire est compilé avec `-march=native` : **recompilez-le sur la machine
qui l'exécute** (`make clean && make`).

Vérifiez l'installation — l'auto-test simule un mini-run de 12 échantillons,
injecte 4 événements connus, les retrouve, **puis compare le lecteur C++ au
lecteur Python** :

```bash
scripts/verifier_installation.sh
```

Vous devez lire `Auto-test : 4/4 evenements retrouves` et
`lecteur C++ vs Python sur 4 echantillons : IDENTIQUES`.

Sans compilation, SpliceFast fonctionne quand même (`--engine auto` retombe
sur le lecteur Python, avec un avertissement), mais beaucoup plus lentement.

**Il vous faut en plus une annotation GTF** correspondant au génome utilisé
pour l'alignement (GENCODE ou Ensembl, `.gtf` ou `.gtf.gz`). Elle sert à
nommer les gènes, à classer les événements et à définir les sites où compter
les reads « bavants » (rétention d'intron). Les noms de chromosomes de
**chaque** BAM sont vérifiés contre ceux du projet : un BAM aligné sur un autre
génome (`21` au lieu de `chr21`) est écarté avec un message clair au lieu de
produire des comptages décalés.

---

## 3. Utilisation sur vos deux runs

### Les scripts prêts à l'emploi

Éditez **une fois** `scripts/config.sh` (dossier de travail, GTF, runs,
panels, nombre de cœurs), puis :

| Script | Rôle |
|---|---|
| `scripts/lancer_splicefast.sh` | estimation de durée puis analyse complète, rétention comprise, en une lecture des BAM |
| `scripts/couverture_introns.sh` | module rétention seul + rapport (à ne lancer que si l'analyse a été faite sans `--retention`) |
| `scripts/inspecter_gene.sh SMC1A [echantillon]` | tableau couverture d'un gène, tous échantillons |
| `scripts/refaire_rapport.sh [padj] [delta]` | rapport seul, avec d'autres seuils |
| `scripts/verifier_installation.sh` | dépendances + auto-test |
| `scripts/diagnostic_disque.sh [minutes]` | test d'endurance du disque (§9) |

Tout ce qui est volumineux (comptages, matrices, caches) va dans le dossier
`--workdir` : placez-le sur un disque qui a de la place (compter ~20 Go pour
150 échantillons). SpliceFast y redirige aussi `TMPDIR`.

### Estimer la durée avant de lancer

```bash
python3 -m splicefast bench \
  --run "/data/BAM et Analyses GenoDiag/RNA-seq_total/novaseq_run4" \
        "/data/BAM et Analyses GenoDiag/RNA-seq_total/novaseq_run5" \
  -w /data/splicefast -t 8
```

Mesure le débit réel du lecteur C++ (Mo/s, à 1 thread et à `-t` threads) sur
les 3 premiers millions de reads, et en déduit la durée totale.

### Lancer l'analyse complète

```bash
python3 -m splicefast all \
  --run "/data/BAM et Analyses GenoDiag/RNA-seq_total/novaseq_run4" \
        "/data/BAM et Analyses GenoDiag/RNA-seq_total/novaseq_run5" \
  --gtf /data/annotation/gencode.v50.basic.annotation.gtf.gz \
  -w /data/splicefast -t 8 -j 1 --retention
```

- `-t` : cœurs utilisés au total ;
- `-j` : BAM lus **simultanément**, les `-t` cœurs étant répartis entre eux
  (`-t 8 -j 2` = 2 BAM à 4 threads de décompression chacun). **`-j 1` sur
  disque externe ou USB** (deux lectures parallèles font chuter le débit d'un
  disque à plateaux), `-j 2` à `4` sur SSD/NVMe ;
- `--retention` : mesure la couverture des introns **pendant la même lecture**
  que le comptage, puis l'intègre au rapport. Sans cette option, le module
  rétention demandera une seconde lecture complète des BAM.

### Ou étape par étape (recommandé la première fois)

```bash
W=/data/splicefast
python3 -m splicefast discover --run /data/.../novaseq_run4 /data/.../novaseq_run5 --out $W/samples.tsv
python3 -m splicefast count    -w $W --samples $W/samples.tsv --gtf /data/annotation/gencode.gtf.gz -t 8 --retention
python3 -m splicefast build    -w $W
python3 -m splicefast analyze  -w $W
python3 -m splicefast retention -w $W        # assemble la couverture mesuree pendant count
python3 -m splicefast report   -w $W
```

`discover` écrit `samples.tsv` (une ligne `nom<TAB>chemin_bam`) : **relisez-le
et corrigez les noms d'échantillons si besoin** avant de lancer `count`.

**Toutes les étapes sont reprenables.** Si la machine s'éteint ou si une étape
échoue, relancez exactement la même commande : les échantillons déjà comptés
**et déjà mesurés pour la couverture** sont ignorés, et `all` saute `build` et
`analyze` s'ils sont déjà terminés. Ajoutez `--force` pour tout refaire.

Si un échantillon n'a pas pu être compté (disque débranché, BAM tronqué,
autre génome), `all` **s'arrête avant `build`** : continuer fabriquerait des
lignes vides dans les matrices et fausserait la référence de toute la cohorte.
Le message indique quoi faire — rebrancher et relancer, ou retirer le BAM
abîmé de `samples.tsv` et relancer avec `--samples`.

### Ajouter des covariables (fortement conseillé)

Le fait que vos échantillons viennent de deux runs différents crée un effet de
lot. Les facteurs latents l'absorbent en général tout seuls, mais vous pouvez
le déclarer explicitement :

```
sample	run	sexe	RIN
PATIENT01	run4	F	7.8
PATIENT02	run5	M	6.2
```

puis `--covariates /data/splicefast/covariables.tsv`. Les colonnes numériques sont
centrées-réduites, les colonnes texte transformées en indicatrices.

---

## 4. Les quatre métriques

Pour chaque jonction *j* on note K = nombre de reads coupés la traversant.

| Métrique | K | N (dénominateur) | Ce qu'elle détecte |
|---|---|---|---|
| **jaccard** | reads de la jonction | reads partageant le donneur **ou** l'accepteur, + reads non coupés au site, − doublons | métrique principale de FRASER 2.0, robuste, peu de faux positifs |
| **psi5** | reads de la jonction | tous les reads partageant le même **donneur** | choix d'accepteur : saut d'exon, accepteur cryptique |
| **psi3** | reads de la jonction | tous les reads partageant le même **accepteur** | choix de donneur : donneur cryptique, pseudo-exon |
| **theta** | reads **coupés** au site | reads coupés + reads **non coupés** chevauchant le site | **efficacité d'épissage** |

**Rétention d'intron = chute de theta.** C'est exactement ce que vous
décriviez : le rapport entre les reads qui « sautent » proprement l'intron et
ceux qui bavent à cheval sur la frontière exon–intron. Un read est compté comme
non coupé s'il couvre le site avec au moins `--anchor` bases (5 par défaut) de
chaque côté, donc un simple chevauchement de bout de read ne compte pas.

> Attention : dans le code R de FRASER, θ est défini dans l'autre sens. Ici la
> convention est « **theta bas = intron retenu** », plus intuitive, et les
> événements sont étiquetés en conséquence (`intron_retention` /
> `efficiency_gain`).

---

## 5. Le modèle statistique

1. **Filtres d'expression** (identiques à FRASER) : au moins un échantillon
   avec K ≥ 10, et 75 % des échantillons avec N ≥ 10.
2. **Facteurs latents** : SVD de logit(ψ) sur un sous-échantillon de jonctions
   bien couvertes → matrice Z (échantillons × q). Ils capturent le batch, le
   sexe, la composition cellulaire, la dégradation de l'ARN… Le nombre q est
   choisi automatiquement (règle du coude, plafonné à S/6).
3. **Valeur attendue** : pour chaque jonction, régression pondérée robuste de
   logit(ψ) sur [1, Z, covariables]. Trois protections contre le défaut
   classique « l'échantillon aberrant apprend son propre modèle » (qui est la
   faiblesse connue de l'auto-encodeur de FRASER) :
   - initialisation sur la **médiane** de la ligne, pas sur un ajustement
     non robuste ;
   - poids de **Tukey redescendants** : au-delà de 3 MAD, le point est *exclu*
     (poids exactement nul), pas seulement atténué ;
   - correction **leave-one-out** exacte par le levier (résidu PRESS) : la
     valeur attendue d'un échantillon est celle qu'on obtiendrait sans lui.
4. **Dispersion** : ρ bêta-binomial par jonction, par la méthode des moments
   puis raffinement par vraisemblance, en **excluant les points aberrants**
   (sinon, avec peu d'échantillons, l'outlier gonfle la dispersion et se cache
   lui-même).
5. **Test bêta-binomial bilatéral**, calculé **en espace log**. C'est
   indispensable : `scipy.stats.betabinom` perd toute précision au-delà de
   p ≈ 1e-11 (vérifié contre `mpmath` à 60 chiffres — sur un cas réel, scipy
   renvoie 1e-11 là où la vraie valeur est 1e-280). Un pré-filtre par
   approximation normale évite le calcul exact là où il ne sert à rien
   (`--screen`, réglé pour ne jamais manquer une p-value < 1e-3).
6. **FDR** : Benjamini-Hochberg **par échantillon** (comme FRASER), puis
   agrégation par gène et second BH entre gènes.

---

### Analyser un échantillon sans polluer la cohorte

Un échantillon peut être **testé sans servir de référence**. Il reçoit ses
p-values comme les autres, mais il ne participe ni aux facteurs latents, ni à
la valeur attendue, ni à la dispersion, ni aux médianes de cohorte. Ses
coordonnées dans l'espace des facteurs sont obtenues par *projection* sur les
axes calculés depuis la référence.

C'est ce qu'il faut pour un séquençage trop peu profond qui sera repassé : sa
couverture intronique erratique gonfle la dispersion estimée (donc fait perdre
de la significativité à toute la cohorte) et relève la rétention « attendue »
(donc masque les vraies rétentions).

```bash
# seuil automatique sur le nombre de reads comptés
python3 -m splicefast reference -w /data/splicefast --min-reads 60e6 --list

# ou une liste explicite
python3 -m splicefast reference -w /data/splicefast \
    --test-only ECH10-PAX ECH11-PAX

# seuls les PAXgene d'au moins 60 M de reads definissent la norme (2.5)
python3 -m splicefast reference -w /data/splicefast --reference-match PAX --min-reads 60e6 --list

# seuls ces échantillons définissent la norme, tous les autres sont testés
python3 -m splicefast reference -w /data/splicefast --reference-file cohorte_propre.txt
```

Les mêmes options existent sur `build` et sur `all`. La commande `reference`
ne recompte rien : elle réécrit `store/reference.json`, puis il suffit de
relancer `analyze` et `report` (quelques minutes), sans relire un seul BAM.

Le nombre de reads est celui compté à l'étape `count` (un read = un
enregistrement, donc ~2 par paire) ; `--list` l'affiche pour chaque
échantillon, ce qui permet de choisir le seuil en connaissance de cause.

**Ce que ça garantit.** Vérifié sur une cohorte simulée : pour les
échantillons de référence, les p-values des quatre métriques et du module de
couverture sont **identiques**, que les échantillons hors référence soient
marqués comme tels ou totalement absents du dossier de travail. Seule la
`padj` bouge de ~1 %, parce que la correction FDR d'un échantillon porte sur
le nombre de jonctions testées, et qu'un échantillon hors référence peut
ajouter ses propres jonctions nouvelles à la matrice. Ces jonctions-là restent
testées : une jonction nouvelle propre à un patient doit l'être, même si ce
patient ne sert pas de référence.

Dans les sorties, la colonne `dans_reference` (oui/non) indique le statut de
chaque ligne, et `n_ech_ref_couverts` compte les échantillons **de référence**
couverts. Le rapport HTML a le filtre correspondant.

**Combien en retirer ?** Chaque échantillon retiré de la référence rétrécit la
cohorte : en dessous de ~30 échantillons de référence, la norme devient
fragile et SpliceFast vous avertit. Avec 147 PAXgene, en retirer 20 ou 30 ne
pose aucun problème.

---

## 6. Ce que vous obtenez — `/data/splicefast/results/`

| Fichier | Contenu |
|---|---|
| **`SpliceFast_viewer.html`** | **visualiseur** : ouvrir, choisir le dossier `results`, cocher des échantillons |
| `<sample>.tsv` | tous les événements d'un échantillon, y compris le contexte des mini-graphiques (JSON) |
| `summary.tsv` | un échantillon par ligne : run, référence, reads, **indice de bruit**, nombre d'anomalies |
| `anomalies.tsv` | une ligne par anomalie (signaux regroupés) |
| `hits_all.tsv` | tous les événements de tous les échantillons (sans le JSON) |
| `panels.tsv` | liste des panels PanelApp (pour le filtre du visualiseur) |
| `results.xlsx` | anomalies, événements, résumé, QC, paramètres (`--no-excel` pour s'en passer) |
| `<sample>.bed` | à glisser dans IGV : chaque événement coloré par type |
| `run_info.json` | version, seuils, date, source PanelApp |

Colonnes principales des TSV, dans l'ordre (PanelApp en tête, phénotypes en dernier) :

| Colonne | Signification |
|---|---|
| `panels`, `panels_green`, `panelapp_best_confidence`, `moi` | PanelApp : panels, panels verts, meilleure confiance, mode de transmission |
| `coordinates` | `chr5:37026328-37027358` : 1re et dernière base de l'intron (1-based, comme FRASER) ; un site seul pour theta |
| `anomaly_type`, `event_type`, `event_detail` | nature de l'anomalie (regroupée) et de chaque signal, ex. `exon extension: +157 bp into intron …` |
| `metrics` | jaccard, psi5, psi3, theta, coverage, combined |
| `split_reads` / `total_reads` | K et N chez le patient |
| `psi`, `psi_expected`, `delta_psi`, `psi_cohort_median` | observé, attendu par le modèle, écart, médiane de la référence |
| `p_value`, `padj_sample`, `padj_gene`, `selected_by` | test, corrections, porte d'entrée (FDR / gene / candidate / concordance) |
| `anomaly_id`, `anomaly_p`, `n_evidence`, `evidence` | regroupement des signaux concordants |
| `retention_p_combined`, `retention_detail` | p de Brown sur les 4 familles de preuves de rétention |
| `n_carriers`, `other_carriers`, `recurrent` | échantillons porteurs du même événement |
| `hypervariable_locus` | HLA / immunoglobulines / récepteurs T |
| `in_reference` | l'échantillon sert-il de norme ? |
| `phenotypes` | phénotypes PanelApp (texte long, laissé à droite) |

### Le visualiseur `SpliceFast_viewer.html`

Une page autonome : pas d'internet, rien n'est envoyé, elle fonctionne dans Edge,
Chrome ou Firefox.

1. Ouvrez la page, cliquez **Open results folder** et choisissez le dossier
   `results` (ou glissez-le sur la page).
2. Cochez à gauche un ou plusieurs échantillons. Pour chacun, la liste indique son run,
   ses reads, son nombre d'anomalies et un badge de **bruit** (orange > 2,
   rouge > 4 fois la médiane de la cohorte).
3. Le tableau affiche **une ligne par anomalie** (décocher pour voir chaque signal), avec
   des lignes compactes. Le locus `chr5:37026328-37027358` tient dans une seule colonne,
   avec un bouton copier. Un clic sur un en-tête trie la colonne.
4. Les filtres :
   - **nom de panel** (liste cherchable de tous les panels PanelApp), avec le niveau de
     confiance (vert, vert + ambre, tout) ;
   - type d'anomalie ;
   - p, |Δψ| et reads minimaux ;
   - **masquer les événements récurrents** (≥ 4 porteurs par défaut) ;
   - masquer les loci HLA/IG/TR ;
   - référence oui/non ;
   - texte libre (gène, phénotype, coordonnées).
5. Un **clic** sur une ligne (puis ↑/↓ pour naviguer) ouvre le détail :
   - la structure du gène avec les arcs de jonction (nombre de reads), les sites theta
     et l'intron couvert ;
   - la **distribution de ψ dans la cohorte**, avec la valeur du patient et l'attendu ;
   - tous les signaux de l'anomalie et la p de rétention combinée ;
   - les autres porteurs ;
   - les panels, avec confiance, MOI et phénotypes ;
   - deux boutons « copier le locus » et « copier la vue IGV » (avec marges).
6. **Sélection et liste de gènes.** La case en début de ligne coche un événement ;
   **Maj + clic** sur une autre case coche toute la plage entre les deux (ou la décoche,
   si la première case était décochée). On peut aussi utiliser Ctrl + clic sur la
   ligne, `Espace` sur la ligne courante, `Maj + ↑/↓`, et la case de l'en-tête pour
   tout cocher. **export gene list** télécharge `splicefast_genes.txt` : les symboles
   des gènes cochés, un par ligne, sans doublon, dans l'ordre du tableau. Si rien
   n'est coché, ce sont les gènes de toutes les lignes affichées.
7. **Export** télécharge en TSV les lignes cochées, ou à défaut toutes les lignes affichées.

### Module dédié : rétention d'intron par la couverture

Les métriques de jonction ne voient la rétention qu'à travers les reads qui
chevauchent la frontière exon/intron — quelques dizaines de reads. Un intron
retenu est en réalité couvert sur toute sa longueur. Ce module mesure donc la
**profondeur** et la normalise très localement :

```
        profondeur moyenne du corps de l'intron
ratio = ---------------------------------------------------------
        moyenne des DEUX bases qui le bordent
        (dernière base de l'exon amont, première base de l'exon aval)
```

On ne prend surtout pas la moyenne de l'exon entier : une bosse de couverture
au milieu d'un exon (3'UTR, exon alternatif, artefact d'amorçage) écraserait
le signal. Les deux bases bordantes sont la référence la plus locale possible,
et le rapport s'interprète directement comme la **fraction de transcrits qui
retiennent l'intron** — à la base bordante, la profondeur compte à la fois les
transcrits épissés et les transcrits retenus.

Les positions exoniques d'**autres** transcrits ou de gènes imbriqués sont
exclues du corps de l'intron : sinon un exon alternatif interne se lirait
comme une rétention.

```bash
python3 -m splicefast retention -w /data/splicefast -t 8
```

Par défaut : **tous les introns annotés**, plus besoin de se limiter à un
panel. Avec `--genes SMC1A SETX` on peut restreindre pour un test rapide.

**La colonne `psi`** du rapport vaut `profondeur_intron / (profondeur_intron +
profondeur_bordures)`, pour rester une proportion testable :

| `psi` | rétention correspondante |
|---|---|
| 0,05 | ~5 % (bruit de fond habituel) |
| 0,10 | ~11 % |
| 0,20 | ~25 % |
| 0,33 | ~50 % |
| 0,50 | 100 % (intron aussi couvert que ses bords) |

**Qui sert de référence.** Deux règles, apprises d'un faux classement réel :

- un échantillon dont la bordure est couverte à moins de `--cov-min-border`
  (10 reads) n'est **ni testé, ni classé, ni utilisé pour ajuster le modèle**.
  Un rapport calculé sur 1 read de bordure n'est pas un résultat : il produit
  des ratios de 10⁹ qui monopolisent les premiers rangs et faussent la norme ;
- un intron est testable dès qu'au moins `--cov-min-samples` échantillons
  (25) y sont correctement couverts — et **non** quand la médiane de toute la
  cohorte est élevée. Un gène à expression très variable (SMC1A : parfois
  5 reads, parfois 200) était sinon écarté d'emblée, alors que 40 échantillons
  suffisent largement à définir la norme.

La référence bordante est mesurée de deux façons, sans surcoût :
la **base unique** qui borde l'intron, et la moyenne des `--border-width`
(25) bases exoniques adjacentes. Le test utilise la seconde — aussi locale
vis-à-vis d'une bosse au milieu de l'exon, mais 25 fois moins bruitée ; la
première reste affichée dans le tableau par gène.

**FDR strict ou liste de relecture.** Avec 60 000 introns testés par
échantillon, il faut `p ≤ 8e-7` pour atteindre `padj ≤ 0,05`. Une rétention
réelle mais modérée — 30 % chez un patient contre 7 % dans la cohorte,
`p = 5e-5` — n'y arrive jamais. Ce n'est pas un défaut du test : c'est le prix
de 60 000 hypothèses, et FRASER se heurte exactement au même mur
(`padjustGene = 1` sur ce genre de ligne).

Le rapport inclut donc, en plus des événements qui passent le FDR, les
**`--cov-top` candidats les plus significatifs de chaque échantillon** (20 par
défaut, `0` pour s'en tenir au FDR strict), pourvu qu'ils passent les filtres
de taille d'effet. Deux colonnes permettent de les lire :

- `padj_echantillon` : FDR sur tous les introns testés. Sévère, et c'est
  voulu — un résultat qui y survit est solide sans discussion ;
- `padj_gene` : correction **à l'intérieur du gène** seulement (Bonferroni sur
  ses introns). C'est la lecture pertinente pour un gène candidat, celle qui
  fait ressortir un événement à `p = 5e-5` avec `padj_gene = 1e-3`. Elle ne
  contrôle **pas** le risque à l'échelle du génome : elle sert à hiérarchiser
  une relecture, pas à conclure seule.

Autres seuils réglables à l'étape `report`, qui ne refait aucune mesure :
`--cov-min-ratio` (0,1 par défaut).

**Comment c'est calculé si vite.** PanDepth alloue un tableau de profondeur
par chromosome et incrémente chaque base de chaque read — imbattable en C, mais
1 Go de mémoire pour chr1. `sfscan` calcule exactement les mêmes nombres sans
matérialiser la profondeur :

- les bornes d'introns, d'exons et des fenêtres bordantes découpent le génome
  en ~3 M de **segments** (les bases bordantes sont isolées dans des segments
  d'une base) ;
- chaque bloc aligné `[a, b)` ajoute sa longueur partielle aux segments de `a`
  et de `b`, et +1/−1 dans un tableau de différences pour les segments
  entièrement couverts. Retrouver ces deux segments demande une recherche
  **par galop** à partir de la position du read précédent : le BAM étant trié,
  la réponse est presque toujours à une ou deux cases, en mémoire chaude ;
- une somme cumulée à la fin donne les bases par segment, donc la profondeur
  moyenne de n'importe quel intervalle par sommes préfixes.

Le même balayage compte, dans la même boucle, les jonctions (table de hachage
à adressage ouvert) et les reads non coupés à chaque site d'épissage : **une
seule lecture du BAM** pour `count` et `retention`. La décompression BGZF, qui
représente l'essentiel du temps restant, est répartie sur les threads par
htslib.

**Validation.** Sur des BAM simulés à partir de GENCODE v50 (chr21 et chr22,
17 516 introns, CIGAR avec soft-clips, insertions, délétions, jonctions
multiples, opérations `=`/`X`, doublons, secondaires, MAPQ faibles), le lecteur
C++ produit des jonctions, des reads non coupés et des profondeurs
**identiques bit à bit** à l'ancien lecteur Python, et 4 000 introns tirés au
hasard ont été vérifiés sans écart contre une profondeur calculée base par
base. Sur un cœur, le débit passe de ~25 Mo/s (Python, deux passes) à
~90 Mo/s (C++, une passe, zlib) et ~180 Mo/s avec libdeflate ; avec plusieurs
cœurs, c'est le disque qui devient la limite. Chaque échantillon mesuré est
mis en cache dans `covcounts/` : un run interrompu reprend où il s'était
arrêté.

Les matrices de profondeur sont stockées en memmap sur le disque
(`store/cov.D.dat`, `store/cov.Bd.dat`) : 147 × 1,2 M introns font 700 Mo par
matrice, qui ne passent jamais par la RAM.

### Instruire un dossier sur un gène précis

Quand on sait quel gène regarder, la correction pour tests multiples n'a plus
lieu d'être : ce qui compte est la position de l'échantillon dans la cohorte.

```bash
python3 -m splicefast gene -w /data/splicefast --gene SMC1A --sample ECH09-PAX_S45_sorted
```

Affiche, intron par intron, la profondeur intronique et bordante de **tous**
les échantillons, le ratio, le z-score robuste, le rang de l'échantillon
demandé parmi les seuls échantillons exploitables, et surtout une **p-value
corrigée sur les seuls introns de ce gène**. C'est le point important : quand
on sait quel gène regarder, la correction multiple ne porte plus que sur une
vingtaine d'introns au lieu de 60 000, et un événement invisible à l'échelle
du génome devient parfaitement documentable. Le tableau complet part dans
`results/gene_SMC1A.tsv`. C'est la commande à utiliser pour documenter un cas
déjà rendu à la famille.

### Annotation PanelApp (tous les panels)

```bash
python3 -m splicefast panelapp --update -w /data/splicefast    # une fois, ~5 min
python3 -m splicefast report   -w /data/splicefast
```

La commande interroge `https://panelapp-aus.org/api/v1/genes/` page par page et écrit
`<workdir>/panelapp/panelapp_genes.tsv.gz`, que `report` utilise automatiquement.
Chaque ligne reçoit :

| Colonne | Contenu |
|---|---|
| `panels` | tous les panels du gène : `Mendeliome [G]; Intellectual disability [A]; …` |
| `panels_green` | les panels où le gène est vert |
| `panelapp_best_confidence` | `green` / `amber` / `red` |
| `moi` | AD, AR, XLR… fusionné depuis les panels verts (à défaut ambres) |
| `phenotypes` | phénotypes dédoublonnés de ces panels |
| `panel_details` | détail par panel (JSON, lu par le visualiseur) |

Confiance PanelApp : 3 et 4 = vert, 2 = ambre, 1 = rouge. Le gène est retrouvé
par son symbole, sinon par son identifiant Ensembl.

**Réseau hospitalier** : `--insecure` si le proxy intercepte le SSL ; ou
télécharger les pages avec `scripts/panelapp_download.ps1` (PowerShell, proxy
Windows) puis `panelapp --from-json DOSSIER`. Les anciens exports TSV
(`--panel Mendeliome.tsv …`) servent encore de secours quand il n'y a pas de cache.

### Un événement connu n'apparaît pas : `why`

Quand un événement repéré ailleurs (FRASER, RT-PCR, un autre run) manque au
rapport, commencez par **`probe`**, qui relit la région dans le BAM :

```bash
python3 -m splicefast probe -w /data/splicefast \
    --sample ECH02-PAX_S61 --region chr5:37025000-37028000
```

Pour chaque jonction et chaque site annoté de la région, `probe` affiche les reads du BAM
ventilés par catégorie (retenus, doublons PCR, secondaires, supplémentaires, MAPQ), le
nombre stocké dans le fichier de comptage, la valeur dans la matrice et la médiane de la
cohorte. Elle conclut ensuite : problème au **comptage** (BAM changé, doublons…), store à
**reconstruire**, ou rien à signaler. Dans ce dernier cas, `why` dit **à quelle étape du
modèle** l'événement a disparu :

```bash
python3 -m splicefast why -w /data/splicefast \
    --sample ECH06-PAX --gene DYNC1H1
# ou par position :
python3 -m splicefast why -w /data/splicefast \
    --sample ECH06-PAX --chrom chr14 --start 102000145 --end 102000285
```

Pour chaque jonction du gène et chaque métrique : k, n, psi observé, psi
attendu, delta, p, padj, et un verdict explicite — `rapporté (FDR)`,
`|delta| = 0.04 < --delta 0.1`, `padj = 0.3 > 0.05 et 87e p-value de
l'échantillon`, ou `jonction non testée : K max cohorte 6 (>= 10 requis)`. La
couverture intronique du gène pour cet échantillon est affichée ensuite.

Les trois cas les plus fréquents, et quoi faire :

| Verdict | Cause | Correctif |
|---|---|---|
| l'échantillon est absent de la liste | il n'a jamais été compté (BAM écarté) | `count` puis `build` |
| `jonction non testée` | trop peu de reads dans la cohorte de référence | `--min-k-one` / `--min-n-quantile` plus bas |
| `padj > seuil` | vrai signal, mais noyé par la correction génome entier | c'est ce que rattrapent les portes `gène` et `candidat` ci-dessous |

### Anomalies : réunir les signaux d'un même événement (2.4)

Un même événement laisse plusieurs traces : une rétention d'intron fait
chuter theta au donneur **et** à l'accepteur, baisser le jaccard de la
jonction canonique et monter la couverture de l'intron. Le rapport les
regroupe en **anomalies** (colonnes `anomalie`, `p_anomalie`, `n_preuves`,
`preuves`) et **trie chaque échantillon par `p_anomalie`**. Le fichier
`anomalies.tsv` (et la première feuille du classeur Excel) donne une ligne par
anomalie.

Pour une rétention, les quatre signaux sont combinés en une p unique
(`p_retention_combinee`, méthode de Brown) : Fisher corrigé de la dépendance
entre métriques, dépendance **mesurée sur la cohorte** (les mêmes reads non
coupés entrent dans theta et jaccard : corrélation ~0,7 à 0,9). Sans cette
correction, Fisher serait 58 fois trop optimiste à p = 1e-4 ; avec elle, le
test reste légèrement conservateur (vérifié sur 20 millions de tirages nuls).

- signaux concordants → p combinée très faible, l'anomalie remonte en tête ;
- couverture isolée, contredite par theta et jaccard (antisens, gène imbriqué)
  → pondérée à la baisse, sans être supprimée ;
- signaux modérés mais concordants → peuvent désormais être rapportés
  (`retenu_par = concordance`), alors qu'aucun ne passait seul ;
- une jonction nouvelle dans le même groupe → type
  `novel_donor+retention_partielle` etc. : c'est le site cryptique qui
  explique la « rétention ».

`psi5` et `psi3` n'entrent pas dans la combinaison : ils ne comparent que des
reads coupés entre eux et sont **aveugles à une rétention pure** (le transcrit
retenu ne produit aucun read coupé ; ψ reste à 1 chez tout le monde).
`--no-harmonise` désactive l'étape.

### Trois portes d'entrée dans le rapport

La colonne `retenu_par` indique pourquoi chaque ligne a été gardée :

- **FDR** : `padj` de l'échantillon ≤ seuil, sur le génome entier ;
- **gène** : à la FRASER2 (`padjustGene`). La p du gène est la meilleure p
  parmi ses jonctions **de taille d'effet suffisante** (|delta| ≥ `--delta`),
  corrigée par le nombre de jonctions testées du gène, puis BH entre gènes.
  **Une seule jonction par gène** entre par cette porte : la plus probante ;
- **candidat** : parmi les `--top-per-sample` meilleures p-values restantes
  (20 par défaut), avec `p ≤ --nominal-p` (1e-4).

Le plafond `--top` (300 lignes par échantillon et par métrique) s'applique
par ordre de priorité — FDR, puis gène, puis candidat — puis de p croissante,
et chaque ligne écartée est signalée dans le journal. `why` utilise
exactement la même fonction de sélection que le rapport.

`--top-per-sample 0` revient au FDR strict. À l'inverse, pour une relecture
exhaustive sur un gène candidat, `why` ignore tous les seuils.

---

## 7. Réglages utiles

| Situation | Option |
|---|---|
| Disque externe lent ou instable | `-j 1` (défaut), `-t 2` à `4` (voir §9) |
| SSD / NVMe, beaucoup de cœurs | `-j 2` à `4` |
| Coupures fréquentes | `--io-retries 30 --io-wait 5` |
| Pas de compilateur | `--engine python` (lent, résultats identiques) |
| Imposer le lecteur C++ | `--engine cpp` (erreur explicite s'il n'est pas compilé) |
| BAM en CRAM | fonctionne avec le lecteur C++ (`--pattern '*.cram'` pour `discover`) |
| Machine avec peu de RAM | `--block-size 5000` (RAM ÷ 4) |
| Trop de jonctions rares | `--min-total 30 --min-samples 3` |
| Restreindre aux gènes codants | `--biotypes protein_coding` |
| Un seul patient à regarder | `--sample PATIENT07` (étape `report`) |
| Plus (ou moins) d'événements à relire | `--html-rows 30000` (15 000 par défaut ; 2 500 lignes dessinées à la fois, le tri et les filtres portent sur la totalité) |
| Rétention d'intron peu sensible | `--retention` avec `all`/`count`, ou `splicefast retention` |
| Événement attendu absent du rapport | `why --sample X --gene Y`, puis `--top-per-sample 50` ou `--nominal-p 1e-3` |
| Échantillon peu profond / douteux | `reference --min-reads 60e6` ou `--test-only NOM` : testé, mais hors norme |
| Nom de run mal deviné | créer `/data/splicefast/runs.tsv` : `echantillon<TAB>run` |
| Recalculer tous les tests exactement | `--screen 1.0` (plus lent, résultat identique aux 3e décimales) |

L'étape `report` peut être relancée autant de fois qu'on veut avec des seuils
différents : elle ne recalcule aucun test, elle relit les p-values déjà
écrites (quelques secondes).

Le lecteur C++ se lance aussi seul, pour vérifier un BAM ou déboguer :
`splicefast/bin/sfscan --help`.

---

## 8. Ordres de grandeur (80 échantillons, 400 Go)

| Étape | Durée | RAM | Disque |
|---|---|---|---|
| `count --retention` | dominée par la lecture des BAM ; le lecteur C++ n'est plus limitant (~150–200 Mo/s par thread de décompression), c'est le **disque** qui fixe le rythme : ~1 h sur un SSD rapide, 2–3 h sur un disque externe à 100 Mo/s | ~100 Mo par BAM lu | ~2 Go |
| `build` | 10–20 min | 500 Mo | ~1,5 Go |
| `analyze` | ~1 h (4 métriques) | ~150 Mo | ~2 Go |
| `retention` (après `count --retention`) | quelques minutes (assemblage) | ~1 Go | ~0,7 Go/matrice |
| `report` | 2–5 min | selon le nombre d'introns testables | quelques Mo |

Ces durées sont des estimations : lancez `bench` pour la mesure réelle sur
votre machine. Il faut bien lire les 400 Go une fois ; la version 2.0 ne les
lit plus qu'**une seule** fois au lieu de deux.

---

## 9. Si le disque se déconnecte ou si un BAM est abîmé

Un disque USB fatigué qui se déconnecte en cours de lecture fait remonter une
erreur d'entrée/sortie (`EIO`). SpliceFast distingue maintenant deux cas :

1. **Le fichier n'est plus accessible** (disque débranché, point de montage
   disparu) : attente, puis nouvelle lecture de cet échantillon — 5 s, 10, 20,
   40… jusqu'à 2 min, et jusqu'à 30 min d'attente du retour du disque. Les
   échantillons déjà terminés sont gardés.
2. **Le fichier est lisible mais son contenu est abîmé** : une seule relecture
   immédiate, puis l'échantillon est mis de côté **sans attendre**. Un BAM
   tronqué (marqueur de fin BGZF absent : copie incomplète, alignement
   interrompu) est détecté dès l'ouverture.

Dans les deux cas le run continue avec les échantillons suivants, la liste des
échantillons écartés est rappelée à la fin, et `all` refuse de construire les
matrices tant qu'il en manque (§3). Relancer la même commande reprend
uniquement ce qui manque.

Tout est tracé dans `/data/splicefast/logs/splicefast.log`.

**Tester le disque avant de lancer des heures de calcul** :

```bash
scripts/diagnostic_disque.sh 20
# et, dans un autre terminal, pour voir les resets USB en direct :
sudo dmesg -w | grep -iE 'usb|reset|I/O error'
```

**À vérifier côté machine** :

- mise en veille automatique du disque : `sudo hdparm -B 254 -S 0 /dev/sdX`
  (ou l'équivalent dans l'utilitaire *Disques*) ;
- mise en veille USB (« autosuspend ») : pour tester, ajouter
  `usbcore.autosuspend=-1` aux options de démarrage du noyau ;
- brancher **directement à l'arrière de la machine**, pas sur un hub, et
  essayer un autre câble : sur un disque à alimentation externe, un câble
  fatigué ou un bloc d'alimentation faiblissant provoque exactement ce
  symptôme sous charge soutenue ;
- si `dmesg` montre des erreurs de secteurs, lancer `sudo smartctl -a /dev/sdX`
  et **copier les données ailleurs** avant d'aller plus loin.

**Et surtout : un seul BAM à la fois sur un disque à plateaux** (`-j 1`, le
défaut). Deux lectures simultanées multiplient les déplacements de tête : le
débit total baisse au lieu de monter.

---

## 10. Limites connues, à avoir en tête

- Les reads non coupés (donc **theta / rétention d'intron**) ne sont comptés
  qu'aux **sites d'épissage annotés du GTF**. C'est un choix délibéré : cela
  permet de tout faire en une seule lecture des BAM. Les jonctions
  *nouvelles*, elles, sont bien détectées partout (jaccard/ψ5/ψ3).
- Les paires de reads qui se chevauchent comptent deux fois sur la zone
  commune. FRASER fait de même. L'effet s'annule au numérateur et au
  dénominateur de ψ, donc il est négligeable.
- Un événement n'est détectable que si le gène est exprimé dans le tissu
  séquencé, et que si le patient est **le seul** à le porter : deux patients
  avec la même anomalie rare s'entre-masquent partiellement.
- Avec moins de ~30 échantillons, les facteurs latents sont peu fiables ;
  utilisez `--n-factors 2` ou `--n-factors 3`.
- Le **run** de chaque échantillon est déduit du chemin de son BAM (le dossier
  qui ressemble à un run : `novaseq_run4`, `nextseq_...`). Si votre arborescence
  ne s'y prête pas, créez `/data/splicefast/runs.tsv` avec deux colonnes
  `echantillon<TAB>run` : il prend le pas sur la déduction automatique. Le run
  détecté est affiché dans le journal à chaque exécution de `report`, vérifiez-le
  d'un coup d'œil.
- Le module de couverture **ignore le brin** : une transcription antisens
  chevauchant un intron gonfle le signal. Le masquage des exons annotés de
  tout transcrit limite le problème sans l'éliminer. Un événement de
  couverture se vérifie dans IGV avant tout.
- Il ne s'agit pas d'un dispositif de diagnostic : tout événement doit être
  vérifié visuellement dans IGV (les fichiers `.bed` sont là pour ça) puis
  confirmé par une méthode indépendante.
