# Fig. 3 | Trajectory analysis identifies distinct state transitions in CD8 and CD4 T cells

> [Back to Repository Overview (EN)](../README.md) · [返回仓库首页 (中文)](../README.zh-CN.md)

![Fig. 3](Combined/M03.png)

[Download PDF](Combined/M03.pdf) · [Original EMF](Combined/M03.emf) · [Primary Shared Inputs](../00_Global_Data/Shared_Inputs/Primary/README.md)

## Key Findings

Monocle2 pseudotime reconstruction reveals that post-transplant `CD8+` T cells (`n = 26,284`) do not progress along a single linear path toward exhaustion, but instead bifurcate from a naive/Tcm root (`State 1`, enriched in HD) into two distinct post-transplant terminal branches: an effector-rich branch (`State 2`, dominated by `S2_CD8_Teff`) and an exhaustion-prone branch (`State 3`, clustering nearly all `CD8+ Tex` together with `Tem` and `S3_CD8_Teff`) (`Fig. 3A–E`). Integrating clonal expansion, tissue distribution, and state transition via STARTRAC identifies **`S2_CD8_Teff`** as the compartment exhibiting the highest combined tumor-reactive clonal profile (`Fig. 3F`). Parallel analysis of `CD4+` T cells (`n = 33,067`) resolves a root naive state (`State 1`) diverging into a post-transplant effector arm (`State 2`, `S2_CD4_Teff`) and a memory/regulatory branch (`State 3`, `Tscm`/`Tcm`/`Treg`), with lower STARTRAC reactivity scores than `S2_CD8_Teff` (`Fig. 3G–L`).

## Panel Descriptions

- **(A–C)** `CD8+` T cell Monocle2 trajectory colored by annotated subset (`A`), pseudotime (`B`), and inferred cellular state (`S1`, `S2`, `S3`; `C`).
- **(D)** Subset composition across `CD8+` trajectory states `S1–S3`.
- **(E)** Stage origin (`HD`, `Pre-HSCT`, `Post-HSCT`) of `CD8+` subsets within `State 2`.
- **(F)** Three-dimensional STARTRAC summary (clonal expansion, tissue migration, and state transition) for `CD8+` T cell states.
- **(G–I)** Parallel `CD4+` T cell trajectory colored by cell type (`G`), pseudotime (`H`), and inferred state (`I`).
- **(J–K)** `CD4+` subset composition across states (`J`) and stage contributions within `State 2` (`K`).
- **(L)** Three-dimensional STARTRAC clonal profile for `CD4+` subsets and states.

## Panel Directory

[A](a_M03a/README.md) · [B](b_M03b/README.md) · [C](c_M03c/README.md) · [D](d_M03d/README.md) · [E](e_M03e/README.md) · [F](f_M03f/README.md) · [G](g_M03g/README.md) · [H](h_M03h/README.md) · [I](i_M03i/README.md) · [J](j_M03j/README.md) · [K](k_M03k/README.md) · [L](l_M03l/README.md)

## Code, Inputs & Reproduction

- **Saved-table verification (no large RDS required)**:
  ```bash
  python 00_Settings/run_saved_plot.py M03 --output /path/to/output_M03
  ```
  Reproduces `CD8` and `CD4` trajectory embeddings (`M03a, c, g, i`), state compositions (`M03d, j`), and stage breakdowns (`M03e, k; S03a–d`) directly from exported coordinates (`CD8cds_*`, `CD4cds_*`) and exact raw-to-saved state mappings (`*_raw_state_to_saved_annotation.csv`).
- **Original archived scripts**: [01_STARTRAC_analysis.R](01_Code/01_STARTRAC_analysis.R) and archived Monocle notes ([01_monocle_analysis.source.txt](00_Preprocess/01_monocle_analysis.source.txt)).
- **Shared STARTRAC outputs**: `00_Global_Data/Shared_Inputs/Primary/Startrac/` (`CD8out.rds`, `CD8data_wide.rds`, `CD8subT.meta.csv`, `CD4subT.meta.csv`).
