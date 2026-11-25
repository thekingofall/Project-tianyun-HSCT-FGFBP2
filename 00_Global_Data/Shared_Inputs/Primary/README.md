# Primary plotting inputs

PBMC files were exported from `scRNA224.rds`: 175,531 cells and 23 samples. Tcell files were exported from `Tcell.rds`: 72,884 cells and 15 subsets. `*_metadata.csv.gz` retains cell IDs, original sample labels and existing annotations; `*_UMAP.csv.gz` retains the stored coordinates. `*_marker_logexpression.csv.gz` contains selected genes from the RNA data slot, without normalization changes. Functional scores and TCR fields remain in Tcell metadata. These gene subsets do not replace full objects for differential expression or new gene queries.

CD8 and CD4 trajectory files preserve fitted coordinates and graph geometry. The raw Monocle State values differ from the saved manuscript-level labels. The two `*_raw_state_to_saved_annotation.csv` tables record the observed cell-ID joins: CD8 raw 1/2/3 → S1, 7 → S2, 4/5/6 → S3; CD4 raw 2 → S1, 1 → S2, 3 → S3. This is a correspondence recovered from saved objects, not a new state assignment.

`Startrac` retains the original clonal inputs and small fitted outputs. `TFnet.csv` is also the byte-identical network file found in the GSE223844 folder and is retained here once.
