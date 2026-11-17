@echo off
REM ===================================================================
REM  SpliceFast 2.5 - reglages communs a tous les .bat : A EDITER UNE FOIS
REM ===================================================================
REM  Dossier de travail : NOUVEAU dossier (les comptages 1.0 ne sont pas
REM  compatibles). Prevoir ~15 Go pour 2 runs.
set WORK=F:\rnaseq\splicefast_2.5
REM  Annotation : la meme que sur Linux de preference
REM  (gencode.v50.basic.annotation.gtf.gz) pour comparer les resultats.
set GTF=E:\annotation\gencode.v44.annotation.gtf.gz
REM  Runs a analyser (recherche des .bam dans les sous-dossiers)
set RUN_A=F:\rnaseq\novaseq_run2
set RUN_B=F:\rnaseq\novaseq_run4
REM  Coeurs utilises, et BAM lus en meme temps (1 sur disque externe/USB)
set THREADS=4
set JOBS=1
REM  Reference (la norme) : seuls les echantillons dont le NOM contient
REM  REF_MATCH ET ayant au moins MIN_READS reads comptes. Les autres sont
REM  analyses (testes) mais ne participent pas a la norme.
set REF_MATCH=PAX
set MIN_READS=60e6
