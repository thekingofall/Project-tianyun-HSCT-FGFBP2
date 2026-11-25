# Supplementary Figure 6 | Clonal lineage connectivity and expression profiles of immune checkpoints across T cell states

> [Back to Repository Overview (EN)](../../README.md) · [返回仓库首页 (中文)](../../README.zh-CN.md)

![Fig. S6](Combined/S06.png)

[Download PDF](Combined/S06.pdf) · [Original EMF](Combined/S06.emf) · [Primary Shared Inputs](../../00_Global_Data/Shared_Inputs/Primary/README.md)

## Key Findings

Alpha- and beta-chain V(D)J alluvial connectivity confirms minimal transition of dominant `S2_CD8_Teff` clones into terminally exhausted `CD8_Tex`, whereas `S3_CD8_Teff` shares extensive dominant clonotypes (`TRAV12-2/TRAJ4`, `TRBV12-4/TRBJ2-7`) with `CD8_Tex` (`Fig. S6A–D`). Consistently, positivity rates for canonical inhibitory receptors (`LAG3`, `PDCD1`, `HAVCR2`, `CTLA4`, and `TIGIT`) are markedly higher in `S3_CD8_Teff` and `S3_CD8_Tex` than in `S2_CD8_Teff`, supporting the low-exhaustion phenotype of `S2_CD8_Teff` (`Fig. S6E`).

## Panel Descriptions

- **(A–D)** V(D)J alluvial plots depicting dominant clonotype sharing for `TRAV–TRAJ` between `S2_CD8_Teff` and `CD8_Tex` (`A`), `S3_CD8_Teff` and `CD8_Tex` (`B`), `S2_CD8_Teff` and `S3_CD8_Teff` (`C`), and `TRBV–TRBJ` between `S2_CD8_Teff` and `S3_CD8_Teff` (`D`).
- **(E)** Percentage of cells positive for inhibitory receptors (`LAG3`, `TIGIT`, `PDCD1`, `CTLA4`, `HAVCR2`) and `TOX` across T cell subsets and inferred states `S1–S3`.

## Panel Directory

[A](a_S06a/README.md) · [B](b_S06b/README.md) · [C](c_S06c/README.md) · [D](d_S06d/README.md) · [E](e_S06e/README.md)

## Code & Reproduction

- **Saved-table verification (`S06e`)**:
  ```bash
  python 00_Settings/run_saved_plot.py S06e --output /path/to/output_S06e
  ```
- **Archived scripts**: [01_marker_positive_fractions.R](01_Code/01_marker_positive_fractions.R), [03_state_specific_VJ_connectivity.R](../../Main_Figure_05/01_Code/03_state_specific_VJ_connectivity.R).
