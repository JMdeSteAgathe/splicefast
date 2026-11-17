# SpliceFast 2.4 — diagnostic des discordances avec FRASER2 et correctifs

## Ce qu'il faut relancer

Seulement `analyze` puis `report` (et `retention` n'a pas à être refait) :

```bash
python3 -m splicefast analyze  -w $run_dir
python3 -m splicefast report   -w $run_dir --panel Mendeliome.tsv Incidentalome.tsv
scripts/verifier_discordances.sh $run_dir
```

Le dernier script contrôle les cinq événements de référence (TBCK, NIPBL,
TUBGCP2, SMC1A, DYNC1H1).

---

## 1. Événements FDR perdus (TBCK, NIPBL, SMC1A/jaccard, TUBGCP2/psi) — bug de sélection

**Cause.** Dans `report.collect_hits` (2.3), quand un échantillon dépassait
`--top` (300) lignes pour une métrique, la place était réservée d'abord aux
lignes « gene » et « candidat » ; les lignes FDR ne recevaient que
`max(0, 300 − n_gene)` places. Au-delà de 300 lignes « gene », **toutes** les
lignes FDR de la métrique étaient supprimées — les plus significatives.

**Preuve sur vos données.** ECH01 : 0 ligne FDR en jaccard/psi3/psi5,
164 en theta (theta n'a pas de porte « gene », il n'était pas touché), et
1 410 lignes « gene ». Toutes les lignes TBCK portent `padj_gene = 2e-25` :
une jonction de TBCK avait p ≈ 1e-27 — celle de FRASER — et a été jetée.
`why` ne le voyait pas : il ne rejouait pas la troncature et annonçait
« rapporté (FDR) ».

**Test de régression** `tests/test_selection.py`, configuration calquée sur
ECH01 : 2.3 → 1 200 lignes dont 1 156 à p > 0,05, vraie anomalie
perdue ; 2.4 → 1 ligne, la vraie anomalie (FDR).

**Correctif.** Une fonction de sélection unique (`report.gate_sample`),
utilisée par le rapport ET par `why`. Plafond appliqué par priorité (FDR,
gène, candidat) puis par p, avec un avertissement dans le journal.

## 2. Spécificité : la porte « gene »

**Cause.** La p du gène pouvait venir d'une jonction à delta minuscule mais
très couverte (AHNAK : p = 1e-63 avec |delta| < 0,1), puis **toutes** les
jonctions du gène à |delta| ≥ 0,1 entraient, quelle que soit leur p.
Sur ECH01 : 935 lignes sur 1 599 à p > 0,05 (jusqu'à 0,95).

**Correctif, à la FRASER2.** La p du gène ne vient que de jonctions de taille
d'effet suffisante, et **une seule jonction par gène** entre par cette porte.
C'est ce qui rend FRASER2 « propre » (aucun autre événement dans le gène).

## 3. ψ attendu = 0,5 — perte de sensibilité ET faux positifs

**Cause.** Les échantillons à N = 0 (gène non exprimé) valaient
logit(0,5) = 0 dans la régression. Dès qu'ils étaient majoritaires (gènes
exprimés chez 25 à 50 % des échantillons : HLA, IG, gènes à expression
variable), la médiane tombait à 0, les poids de Tukey excluaient tous les
échantillons informatifs et ψ attendu valait 0,5 pour tout le monde.
Sur vos TSV : 77 et 55 lignes avec ψ attendu loin de la médiane de cohorte,
beaucoup à exactement 0,5.

**Démonstration** (`tests/test_modele.py`) : outlier à ψ = 0,6 contre 0,03
dans la cohorte — 2.3 : ψ attendu 0,50, dispersion 0,6, **p = 0,91** ;
2.4 : ψ attendu 0,035, **p = 3e-102**.

**Correctif.** Seuls les échantillons de référence avec N ≥ `--min-cov`
ajustent le modèle (le module couverture le faisait déjà).

## 4. Facteurs latents qui apprennent les outliers

**Cause.** Un échantillon très atypique pouvait porter à lui seul un axe de la
SVD ; la régression « expliquait » alors ses anomalies. Sur l'auto-test, le
jaccard de la jonction canonique de SIM09 (ψ 0,12 contre 0,9) sortait avec
ψ attendu 0,013, delta +0,11 : anomalie absorbée. Dans `tests/test_modele.py`,
le premier axe 2.3 est l'échantillon aberrant (|r| avec le lot = 0,20) au
lieu du lot (2.4 : |r| = 1,00).

**Correctif.** Winsorisation des lignes à 4 MAD avant la SVD, et rejet des
axes portés par moins de max(3, 2 % de la référence) échantillons effectifs.
Chaque axe rejeté est signalé dans le journal avec l'échantillon qui le
porte : **c'est un indicateur de QC** (échantillon dégradé ou très atypique).
Après correctif : SIM09 jaccard delta −0,61, p = 1e-24.

## 5. Harmonisation : anomalies et p combinée de rétention

Nouvelle étape (`harmonize.py`), voir README « Anomalies ». Chaque ligne reçoit
`anomalie`, `p_anomalie`, `n_preuves`, `preuves` ; les échantillons sont triés
par `p_anomalie` ; `anomalies.tsv` et la première feuille Excel donnent une
ligne par anomalie.

Rétention : couverture, theta donneur, theta accepteur et jaccard combinés en
p unilatérales par la méthode de Brown, dépendance mesurée sur la cohorte,
variance gonflée de 10 %. Calibration sur 20 M de tirages nuls corrélés
(ρ jusqu'à 0,9) : Brown ×0,55 à ×0,74 de α entre 1e-3 et 1e-6 (conservateur) ;
Fisher naïf ×58 à 1e-4.

Auto-test : la rétention SIM05 passe de p = 1e-28 (meilleur signal seul) à
4e-32 (4 preuves) ; le donneur cryptique SIM09 est étiqueté
`novel_donor+retention_partielle`.

## 6. Divers

- `result.html` : quota garanti par échantillon (en 2.3, un tri global unique
  sur 15 000 lignes pouvait faire sortir de la page un événement à p = 1e-14
  d'un patient) ; tri par défaut sur `p anomalie` ; colonnes p, anomalie,
  preuves.
- `why` : coordonnées en 1-based (comme le rapport et FRASER) ; verdicts
  « écarté par --top » et « une autre jonction du gène est plus probante ».
- `--no-harmonise` pour retrouver un rapport sans regroupement.

---

## Réponses aux questions

**SMC1A — pourquoi psi5/psi3 ne s'allument pas ?** Ils comparent des reads
coupés entre eux. Un transcrit qui retient l'intron n'en produit aucun : si
aucune autre jonction ne partage le donneur ou l'accepteur, ψ5 = ψ3 = 1 chez
le patient comme dans la cohorte. Aveugles par construction. Jaccard, lui,
intègre les reads non coupés : 116 / (116 + 251 + 140) ≈ 0,23 contre ~0,85
attendu. Il aurait dû sortir : bug n° 1.

**TUBGCP2 — pourquoi seulement la couverture ?** Chez FRASER, la proportion
de reads non coupés (0,23) est **sous** son 99e centile (0,28) : ce n'est pas
un excès de rétention. En revanche ~800 reads coupés (1146 − 85 − 264)
partagent le donneur ou l'accepteur : une jonction concurrente, probablement un
site cryptique intra-intronique (ce qui explique aussi la hausse de
couverture). theta ne peut pas la voir : ces reads comptent comme « épissés »
au site. psi5/psi3/jaccard l'avaient vue — bug n° 1. À vérifier dans IGV ;
en 2.4 ce cas doit sortir en `novel_*+retention_partielle`.

**NIPBL.** Intron annoté : rien ne filtre les introns annotés. Jaccard
55/105 avec excès de reads non coupés : signal FDR très probablement supprimé
par le bug n° 1. Le script de vérification le confirmera.

**DYNC1H1 absent en 2.2, présent en 2.3.** Inexplicable par le code à store
identique (une ligne psi5 à p = 1e-14 sortait forcément en 2.2). Hypothèses :
store ou référence différents entre les deux analyses, ou lecture dans
`result.html`, limité en 2.2 aux 15 000 meilleures lignes de la cohorte — en
2.3 le bug n° 1 a libéré beaucoup de place dans ce top 15 000. Un
`grep DYNC1H1` dans le `hits_all.tsv` de la 2.2 tranchera.

## Limites et ce qui reste à valider sur données réelles

- Validé sur simulations et tests unitaires seulement : les chiffres de
  spécificité sur votre cohorte (lignes par échantillon, lignes à p > 0,05)
  sont à mesurer après relance.
- Le rejet d'axes latents peut écarter un vrai petit lot (< 2 % de la
  référence, soit < 7 échantillons sur 340) : surveiller le journal.
- Les anomalies non-rétention (saut d'exon, pseudo-exon…) sont regroupées mais
  pas combinées statistiquement : leur `p_anomalie` est la meilleure p du
  groupe.
- La combinaison suppose que la dépendance entre métriques est la même pour
  tous les introns ; elle est estimée sur 60 échantillons tirés au hasard.
