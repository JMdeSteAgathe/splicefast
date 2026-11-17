# SpliceFast 2.7

## 1. SYNGAP1 (ECH07) : détecté, mais masqué par deux filtres trop larges

**Constat.** FRASER rapporte une rétention de l'intron chr6:33438564-33438774
(Jaccard 82/433, Δψ = −0,75). Le fichier 2.6 de l'échantillon contient bien
l'événement, et il est significatif après FDR sur trois signaux indépendants :

| signal | k/n | p |
|---|---|---|
| Jaccard de l'intron | 96/152 | 6e-54 |
| couverture intronique (module retention) | — | 5e-37 |
| theta (reads non coupés aux sites) | — | 7e-33 |

Le modèle n'est donc pas en cause : l'anomalie était **masquée dans le
visualiseur**, qui cache par défaut les loci hypervariables et les événements
récurrents. Deux erreurs l'y faisaient tomber.

1. **« Locus hypervariable » défini par région.** SYNGAP1 se trouve dans la région
   du CMH (chr6:28,5-33,5 Mb), et la 2.5/2.6 marquait comme hypervariable *tout*
   gène de cette région. Dans votre cohorte, **197 gènes non immunitaires**
   étaient concernés, soit 17 397 anomalies : SYNGAP1, DAXX, TAP1, COL11A2, TRAF3,
   etc. Désormais, seul le **nom du gène** est pris en compte : HLA-*, gènes
   d'immunoglobulines (IGH/IGK/IGL V-D-J-C), récepteurs T (TRA/TRB/TRG/TRD
   V-D-J-C), KIR, LILR, MICA/MICB.
2. **Récurrence trop permissive.** L'anomalie était comptée chez 4 à 5
   échantillons, dont l'un très bruyant (ECH08, 1 446 anomalies), alors
   que vous n'en aviez retrouvé qu'un, peu convaincant. Un autre échantillon n'est
   plus compté comme porteur que si les trois conditions suivantes sont réunies :
   - p ≤ 1e-4 sur **le même** signal (jonction, site ou intron) ;
   - même sens ;
   - |Δψ| ≥ max(seuil, **la moitié du Δψ de l'échantillon étudié**). Un Δψ de
     −0,2 chez un autre ne fait donc plus de lui un porteur d'un Δψ de −0,75.

   Les échantillons très bruyants (indice de bruit > 4) sont **listés à part et
   ne comptent pas** : nouvelle colonne `n_carriers_noisy`, et mention
   `[+N noisy]` dans `other_carriers`.

**À relancer : `report` seulement** (quelques minutes). Comptages, matrices et
tests restent valables.

```bash
python3 -m splicefast report -w $run_dir
```

Contrôle : ouvrez le visualiseur, cochez ECH07, tapez `SYNGAP1`. Si la
ligne n'apparaît pas, décochez `hide recurrent` : le panneau de droite affiche
alors les porteurs et leurs scores, ce qui permet de juger.

Test de non-régression : `tests/test_recurrence.py`.

## 2. Visualiseur

- **Colonnes réordonnées** pour garder phénotypes et mode de transmission
  visibles quand le panneau de droite est ouvert : Sample, Gene, Phenotypes,
  Inheritance, Anomaly, Signals, Reads, ψ, Δψ, p, padj, Carriers, Panels, puis
  Locus en dernier. La hauteur des lignes reste fixe : les phénotypes sont
  tronqués, et le texte complet s'affiche au survol.
- **Panneau des échantillons repliable** : bouton `«` / `»` en haut à gauche.
  L'état est mémorisé dans le navigateur.
- **Bouton `all`** : il sélectionne exactement les échantillons affichés par le
  filtre, que celui-ci porte sur le nom ou sur le run (taper `run1` puis `all`
  fonctionne désormais).
- **Qualité de l'échantillon** en tête du panneau de droite : indice de bruit
  (global, puis jonctions / couverture), nombre de reads, fraction de reads
  épissés, nombre d'anomalies (dont FDR), référence ou test seul.
- **Tableau des porteurs** : l'échantillon courant (surligné) figure à côté de
  chaque autre porteur, avec Δψ, k/n et p. Les échantillons bruyants sont grisés
  et ceux hors référence sont signalés. Un clic sur un nom ajoute cet échantillon
  à la sélection, et ses lignes apparaissent dans le tableau. La colonne
  `Carriers` affiche aussi `+N` pour les porteurs bruyants non comptés.
- **Schéma du gène** : le nombre de reads est placé au-dessus de l'arc, sur un
  fond masquant, et reste lisible en mode sombre.

## 2.7.1 — largeur des colonnes

- Les colonnes gardent leur largeur. En 2.7.0, le tableau s'élargissait jusqu'au
  plus long texte (`width:max-content`), si bien que la colonne Phenotypes
  repoussait toutes les autres hors de l'écran. Le texte trop long est maintenant
  coupé par « … », et le texte complet s'affiche au survol.
- La colonne Phenotypes affiche environ 64 caractères. Sa largeur est mesurée dans
  la police de votre navigateur.
