# =============================================================================
#  Reglages communs a tous les scripts : A EDITER UNE FOIS.
# =============================================================================
WORK=/data/splicefast                         # dossier de travail (gros disque)
GTF=/data/annotation/gencode.v50.basic.annotation.gtf.gz
PANEL_DIR=$HOME/panelapp_g2p                  # Incidentalome.tsv, Mendeliome.tsv
RUNS=(
  "/data/BAM et Analyses GenoDiag/RNA-seq_total/novaseq_run4"
  "/data/BAM et Analyses GenoDiag/RNA-seq_total/novaseq_run5"
)
THREADS=$(nproc)      # coeurs utilises
JOBS=1                # BAM lus en meme temps : 1 sur disque externe/USB, 2-4 sur SSD
PANELS=("$PANEL_DIR/Incidentalome.tsv" "$PANEL_DIR/Mendeliome.tsv")
