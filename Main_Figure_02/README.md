# Fig. 2 | Remodelling of T cell functional states and the TCR repertoire post-allo-HSCT

> [Back to Repository Overview (EN)](../README.md) · [返回仓库首页 (中文)](../README.zh-CN.md)

![Fig. 2](Combined/M02.png)

[Download PDF](Combined/M02.pdf) · [Original EMF](Combined/M02.emf) · [Primary Shared Inputs](../00_Global_Data/Shared_Inputs/Primary/README.md)

## Key Findings

High-resolution re-clustering of **72,884 T cells** identifies **15 transcriptionally distinct subsets** spanning naive, memory (`Tscm`, `Tcm`, `Tem`), cytotoxic effector (`Teff`), terminally exhausted (`Tex`), regulatory (`Treg`), and innate-like (`Trm-like`, `NKT`, `MAIT`, `γδT`) compartments (`Fig. 2A–B`). Despite numerical recovery by month 6, post-transplant reconstitution exhibits persistent naive T cell depletion and pronounced expansion of cytotoxic `CD8+ Teff` and `CD8+ Tex` populations (`Fig. 2C`). Functionally, reconstituting T cells display elevated activation, cytotoxicity, and terminal exhaustion module scores peaking at `Post-1M` (`Fig. 2D`). Paired scTCR-seq demonstrates sharp clonal contraction and persistently reduced repertoire diversity (`Shannon`, `Inv. Simpson`, `Chao`, `ACE`, and `Inv. Pielou`) throughout the 6-month window (`Fig. 2E–F`).

## Panel Descriptions

- **(A)** UMAP projection of the 15 transcriptionally defined T cell states.
- **(B)** Canonical lineage and functional marker dot plot (`IL7R`, `CCR7`, `TCF7`, `SELL`, `GNLY`, `PRF1`, `GZMH`, `GZMK`, `LAG3`, `TIGIT`, `PDCD1`, and unconventional lineage markers).
- **(C)** Longitudinal composition of T cell subsets across `HD`, `Pre-HSCT`, `Post-1M`, `Post-3M`, and `Post-6M`.
- **(D)** Violin-box plots of activation, cytotoxicity, and terminal exhaustion module scores across stages (two-sided Wilcoxon rank-sum test).
- **(E)** Paired scTCR-seq clonotype abundance curves showing repertoire skewing after allo-HSCT.
- **(F)** Longitudinal TCR repertoire diversity indices (`Shannon`, `Inv. Simpson`, `Chao`, `ACE`, `Inv. Pielou`).

## Panel Directory

[A (`a_M02a`)](a_M02a/README.md) · [B (`b_M02b`)](b_M02b/README.md) · [C (`c_M02c`)](c_M02c/README.md) · [D (`d_M02d`)](d_M02d/README.md) · [E (`e_M02e`)](e_M02e/README.md) · [F (`f_M02f`)](f_M02f/README.md)

## Code, Inputs & Reproduction

- **Saved-table verification (no large RDS required)**:
  ```bash
  python 00_Settings/run_saved_plot.py M02 --output /path/to/output_M02
  ```
  Generates `M02a` (15-subset UMAP), `M02b` (canonical marker dot plot), `M02c` (subset composition), and `M02d`/`S02d` (functional score violins) directly from `00_Global_Data/Shared_Inputs/Primary/Tcell_*.csv.gz`.
- **Original archived scripts**: [01_scRepertoire_analysis.R](00_Preprocess/01_scRepertoire_analysis.R), [01_primary_landscape_and_TCR.R](../Main_Figure_01/01_Code/01_primary_landscape_and_TCR.R), and archived score notes (`01_Code/*.source.txt`).
- **Primary server objects**: `/data5/tiany/HSCT/Result/Tcell.rds` (72,884 cells) and `/data5/tiany/HSCT/Result/final_vdj.rda` (needed for full `scRepertoire` bootstrapping in `M02e–f`).
