# SpliceFast 2.7 sous Windows — installation et premier run

## Mise à jour depuis la 2.6

Remplacez le dossier `splicefast` (gardez votre `config.bat`), puis lancez
`refaire_rapport.bat` : seul le rapport est refait, en quelques minutes (voir
`CHANGELOG_2.7.md`).

Depuis la 2.5 ou une version antérieure, relancez plutôt `lancer_splicefast.bat`.
Rien n'est recompté, mais les matrices et les tests sont refaits (≈ 1 h).

## 1. Installer (5 minutes)

Rien à compiler : `splicefast\bin\sfscan.exe` est fourni.

1. Fermez toute fenêtre SpliceFast ouverte.
2. Dans `%USERPROFILE%\Documents\code\psl\rnaseq`, **supprimez** :
   - le dossier `splicefast` (version 1.0 : d'anciens modules resteraient sinon) ;
   - les anciens `.bat` (`couverture_introns.bat`, `retention_panels.bat`, etc.).

   Ne touchez pas à `rnu12`.
3. Décompressez `splicefast-2_7-windows.zip` **dans ce même dossier**. Vous obtenez :
   ```
   rnaseq\
     splicefast\            (le programme, avec bin\sfscan.exe)
     config.bat             <- seul fichier a editer
     lancer_splicefast.bat
     refaire_rapport.bat
     mise_a_jour_panelapp.bat, panelapp_download.ps1
     sonder_region.bat, pourquoi.bat, inspecter_gene.bat
     ouvrir_resultats.bat, verifier_installation.bat, diagnostic_disque.bat
     README_FR.md, CHANGELOG_2.7.md
     rnu12\                 (inchange)
   ```
4. Si Windows a marqué les fichiers comme téléchargés : clic droit sur le zip,
   **Propriétés**, **Débloquer**, *avant* de décompresser. Sans cela, l'antivirus peut
   bloquer `sfscan.exe`. Si c'est le cas, SpliceFast continue avec le lecteur Python,
   mais beaucoup plus lentement.
5. Double-cliquez sur **`verifier_installation.bat`**. Il faut lire :
   - `sfscan 2.1.0 (htslib 1.21)` ;
   - `Auto-test : 4/4 evenements retrouves` ;
   - `lecteur C++ vs Python sur 4 echantillons : IDENTIQUES`.

   Python, numpy et scipy sont ceux déjà installés pour la 1.0 (Python 3.14).

## 2. Régler `config.bat` (une fois)

```bat
set WORK=F:\rnaseq\splicefast_2.5
set GTF=E:\annotation\gencode.v44.annotation.gtf.gz
set RUN_A=F:\rnaseq\novaseq_run2
set RUN_B=F:\rnaseq\novaseq_run4
set THREADS=4
set JOBS=1
set REF_MATCH=PAX
set MIN_READS=60e6
```

- **WORK** : un *nouveau* dossier. Les comptages de la 1.0 ne sont pas compatibles.
  Prévoyez ~15 Go.
- **GTF** : c'est votre ancien chemin. Pour comparer avec Linux, mieux vaut la même
  annotation (`gencode.v50.basic.annotation.gtf.gz`) copiée sur le PC.
- **THREADS / JOBS** : 4 et 1 conviennent à un disque externe. Sur SSD, `JOBS=2`.
- **REF_MATCH / MIN_READS** : seuls les échantillons dont le nom contient `PAX`
  (sans tenir compte de la casse) **et** qui ont au moins 60 millions de reads
  comptés forment la référence. Les autres sont analysés, mais ne participent
  pas à la norme. Pour vérifier la liste : `voir_reference.bat`, ou la colonne
  `in_reference` de `summary.tsv`.

## 3. Lancer

Double-cliquez sur **`lancer_splicefast.bat`** :
1. il estime la durée et attend une touche ;
2. il télécharge PanelApp (tous les panels) s'il n'est pas encore en cache ;
3. il enchaîne comptage et couverture (une seule lecture des BAM), matrices, tests,
   puis rapport ;
4. à la fin, il ouvre le visualiseur. Cliquez **Open results folder** et choisissez
   `F:\rnaseq\splicefast_2.5\results`.

Après une coupure (disque, veille), relancez le même `.bat` : le travail déjà fait est
gardé. Le journal se trouve dans `F:\rnaseq\splicefast_2.5\logs\splicefast.log`.

## 4. Ensuite

| Besoin | Fichier |
|---|---|
| Changer les seuils (rapport seul, quelques minutes) | `refaire_rapport.bat 0.1 0.05` |
| Rafraîchir PanelApp | `mise_a_jour_panelapp.bat` puis `refaire_rapport.bat` |
| Un événement vu dans IGV est absent | `sonder_region.bat ECH02-PAX_S61 chr5:37025000-37028000` |
| … puis, si le comptage est bon | `pourquoi.bat ECH02-PAX_S61 NIPBL` |
| Couverture d'un gène, tous échantillons | `inspecter_gene.bat SMC1A` |
| Rouvrir les résultats | `ouvrir_resultats.bat` |
| Voir qui fait partie de la référence | `voir_reference.bat` |
| Changer la référence (éditer `config.bat`), sans recompter | `changer_reference.bat` |

**Premier contrôle conseillé** : novaseq_run2 contient ECH02. Dans le
visualiseur, cochez-le et tapez `NIPBL`. L'anomalie attendue est
`exon_extension (partial intron retention)` en chr5:37026328-37027358, avec la jonction
nouvelle chr5:37026485-37027358.

**Ce test ne porte que sur 2 runs.** La cohorte de référence est plus petite que sur
Linux, ce qui donne une norme un peu moins précise. Un événement récurrent sur Linux peut
ne pas l'être ici, et inversement.

## 5. PanelApp derrière le proxy de l'hôpital

`mise_a_jour_panelapp.bat` essaie d'abord en Python, puis via PowerShell, qui utilise le
proxy et les certificats de Windows. Si les deux échouent, l'analyse continue sans
annotation de panels. Vous pourrez relancer `mise_a_jour_panelapp.bat` plus tard, puis
`refaire_rapport.bat`.
