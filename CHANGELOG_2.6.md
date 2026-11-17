# SpliceFast 2.6

## 1. Faux positifs à ψ = 1 (YY1 et consorts) : un bug de matrice, pas de ψ

**Constat (ECH05, YY1).** La jonction chr14:100262467-100262574 est
rapportée avec ψ5 = 10/10. Or le BAM compte 1 142 reads sur la jonction
majeure chr14:100262467-100274697, qui partage le même donneur. ψ5 aurait dû
valoir 10/1 152 ≈ 0,9 %.

**Cause.** La définition de ψ n'est pas en cause : la jonction majeure était à
**0 read dans la matrice** pour cet échantillon (ligne `annotated_intron` de YY1 :
`split_reads 0`, `total_reads 10`). Dans le même fichier, **233 introns annotés**
sont dans ce cas : leur ψ attendu est ≈ 1, mais leur nombre de reads stocké est nul (PTPN6,
FAS, POMP, DHRS7…). Chacun produit une paire de faux positifs : l'intron
canonique « perdu » et la jonction minoritaire voisine « à ψ = 1 ».

Lors de la construction des matrices (`build`), les jonctions sont rangées
**par cluster** de sites partagés. La liste n'est donc plus triée par position,
alors que le remplissage des comptages la parcourait par recherche
dichotomique, ce qui suppose un tri. Quand un cluster contient une jonction
située après le début d'un autre cluster (gènes imbriqués, sites alternatifs
partageant un accepteur), une partie des jonctions n'était pas retrouvée : elles
recevaient 0 read, sans aucun message. Le défaut existait dès la 2.0, mais le
regroupement exact de la 2.5 l'a beaucoup aggravé, car les clusters sont plus
grands et le désordre plus fréquent.

**Correctif.**
- La recherche se fait maintenant sur une copie triée des clés.
- `build` vérifie, pour chaque échantillon, que **tous** les reads des jonctions
  retenues sont bien placés dans la matrice ; à la moindre perte, il s'arrête
  au lieu de produire une matrice fausse.
- Le store porte un numéro de format : `all` reconstruit automatiquement un store
  antérieur à la 2.6, et `analyze` le signale.
- Test de régression : `tests/test_modele.py` [4].

**Faut-il en plus « corriger » ψ ?** Non, une fois les comptages justes. ψ5 et ψ3
sont les métriques de FRASER, et le Jaccard combine déjà les deux extrémités
avec les reads non coupés. Une jonction dont une extrémité n'est utilisée par
personne d'autre a bien ψ = 1 chez tout le monde, et n'est donc jamais aberrante.
Les ψ = 1 aberrants venaient exclusivement des jonctions voisines mises à zéro.
Ajouter une pondération ad hoc aurait masqué le bug au lieu de le corriger, et
coûté de la sensibilité sur de vrais événements minoritaires.

**À relancer : `build`, `analyze` puis `report`, sans recompter.**
`lancer_splicefast.bat` (Windows) ou `all` (Linux) le font tout seuls :
comptages et couverture sont réutilisés, et seules les matrices et les tests sont
refaits (environ 1 h pour 170 échantillons).

## 2. Gain d'efficacité : canonique ou non

`efficiency_gain` est maintenant scindé :
- `efficiency_gain_canonical` : le site appartient à un intron d'un transcrit
  `Ensembl_canonical` du GTF ;
- `efficiency_gain_non-canonical` : dans le cas contraire.

Une nouvelle colonne `ensembl_canonical` (yes/no) est remplie pour toutes les
jonctions et tous les sites. Elle est lue dans le GTF du projet (une fois, avec
mise en cache). Si le GTF n'a aucun tag `Ensembl_canonical`, ce sont les
transcrits `MANE_Select` qui servent de référence.

## 3. Visualiseur

- Un seul bouton, **copy coordinates**.
- Filtre **PanelApp** : `no filter` (par défaut), `in a panel (any confidence)`,
  `green or amber`, `green only`. Le filtre s'applique aux panels choisis, ou à
  tous si aucun n'est choisi : « green only » ne garde que les gènes verts dans au
  moins un panel. Les gènes absents de PanelApp sont écartés dès qu'un niveau est
  choisi.
- Schéma du gène : les exons de **tous** les signaux de l'anomalie sont
  réunis. Un site theta (gain d'efficacité, rétention) est toujours dessiné
  entre ses deux exons, même seul.
