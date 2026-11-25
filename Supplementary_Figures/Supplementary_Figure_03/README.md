# Supplementary Figure 3 | Compositional analysis of T cell developmental states and external validation of transcriptomic signatures

> [Back to Repository Overview (EN)](../../README.md) · [返回仓库首页 (中文)](../../README.zh-CN.md)

![Fig. S3](Combined/S03.png)

[Download PDF](Combined/S03.pdf) · [Original EMF](Combined/S03.emf) · [Primary Shared Inputs](../../00_Global_Data/Shared_Inputs/Primary/README.md)

## Key Findings

Stage-origin breakdowns confirm that root `State 1` in both `CD8+` and `CD4+` trajectories derives predominantly from `HD` samples, whereas terminal `State 3` is enriched in post-transplant (and for `CD4+`, `Pre-HSCT`) cells (`Fig. S3A–D`). External tumor-reactive `NeoTCR` signature scoring further highlights enhanced antigen-responsive effector potential within post-transplant effector compartments (`Fig. S3E–F`).

## Panel Descriptions

- **(A–B)** Stage contributions (`HD`, `Pre-HSCT`, `Post-HSCT`) to `CD8+` T cell subsets within `State 1` (`A`) and `State 3` (`B`).
- **(C–D)** Stage contributions to `CD4+` T cell subsets within `State 1` (`C`) and `State 3` (`D`).
- **(E–F)** Violin plots of external `NeoTCR` signature scores (Lowery et al., *Science* 2022) across `CD8+` (`E`) and `CD4+` (`F`) state-specific subsets.

## Panel Directory

[A](a_S03a/README.md) · [B](b_S03b/README.md) · [C](c_S03c/README.md) · [D](d_S03d/README.md) · [E](e_S03e/README.md) · [F](f_S03f/README.md)

## Code, Inputs & Known Source Boundary

- **Saved-table verification (`S03a–d`)**: `python 00_Settings/run_saved_plot.py M03 --output /path/to/output_M03`.
- **Known source boundary (`S03e–f`)**: The exact `NeoTCR` scoring script and saved score column were not present in the archived `Tcell.rds` metadata or scripts on server 252; the submitted manuscript panels are retained verbatim in `Combined/S03.emf`.
