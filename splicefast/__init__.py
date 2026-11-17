"""
SpliceFast - detection d'evenements d'epissage aberrants (type FRASER)
concue pour de gros runs sur une machine Linux modeste, sans droits root,
sur des runs de 70-80 echantillons / plusieurs centaines de Go de BAM.

Metriques calculees (equivalentes a FRASER / FRASER 2.0) :
  psi5    : usage du site donneur   (k / somme des jonctions partageant le donneur)
  psi3    : usage du site accepteur (k / somme des jonctions partageant l'accepteur)
  theta   : efficacite d'epissage du site = reads coupes / (coupes + non coupes)
            -> une BAISSE de theta = retention d'intron
  jaccard : Intron Jaccard Index (FRASER 2.0), robuste, combine les 3 signaux

Auteur : genere pour un usage diagnostique de recherche.
"""

__version__ = "2.7.1"
__all__ = ["__version__"]
