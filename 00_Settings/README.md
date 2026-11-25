# Execution, Reproducibility & Provenance Notes / 执行、复现与溯源说明

> [Back to Repository Overview (EN)](../README.md) · [返回仓库首页 (中文)](../README.zh-CN.md)

## 1. Panel Index & Locked Palettes

- [00_Panel_Code_Index.tsv](00_Panel_Code_Index.tsv) maps all **71 submitted panels** across 5 Main Figures and 6 Supplementary Figures to their panel folders, captions, archived scripts, and provenance status.
- [GLOBAL_PALETTE.md](GLOBAL_PALETTE.md) and [global_palette.tsv](global_palette.tsv) record the exact label-to-HEX color bindings extracted from the submitted primary and external validation scripts.
- [01_File_Provenance.tsv](01_File_Provenance.tsv) records SHA256 checksums for all 318 archived source assets copied from server 252.

## 2. Verified Saved-Table Reproduction (`run_saved_plot.py`)

`run_saved_plot.py` resolves the repository root from its own location, stages compact input tables into a fresh output directory, sanitizes server-specific paths and non-interactive graphics device settings without modifying archived originals, and executes R end-to-end.

All **11 saved-table reproduction jobs** were executed and verified end-to-end with **R 4.1.0** (`R_sessionInfo.txt`):

| Job ID | Manuscript Scope | Input Source (Compact Repository Tables) | Script |
|---|---|---|---|
| `M01` | Fig. 1B–E (PBMC UMAP, markers, composition, stage density) | `Shared_Inputs/Primary/PBMC_*.csv.gz` | `00_Settings/03_reproduce_primary_panels.R` |
| `M02` | Fig. 2A–D & S2D (T cell UMAP, markers, composition, functional scores) | `Shared_Inputs/Primary/Tcell_*.csv.gz` | `00_Settings/03_reproduce_primary_panels.R` |
| `M03` | Fig. 3A–E, G–K & S3A–D (CD8/CD4 trajectories & state breakdowns) | `Shared_Inputs/Primary/CD8cds_*`, `CD4cds_*` | `00_Settings/03_reproduce_primary_panels.R` |
| `S06e` | Fig. S6E (Inhibitory receptor & TOX positivity across T cell states) | `Shared_Inputs/Primary/Tcell_*.csv.gz` | `00_Settings/03_reproduce_primary_panels.R` |
| `M04h` | Fig. 4H (GSE216187 16-marker CR vs REL heatmap) | `Shared_Inputs/GSE216187_effector_DE/` | `14_heatmap_markers16_compact_v4.R` |
| `M04i` | Fig. 4I (GSE222633 CR vs REL longitudinal trajectories) | `Shared_Inputs/GSE222633_CR_REL/` | `02_redraw_CR_REL.R` |
| `S05b` | Fig. S5B (GSE216187 24-celltype canonical marker heatmap) | `Shared_Inputs/GSE216187_markers/` | `11_render_marker_heatmap_final.R` |
| `S05c` | Fig. S5C (GSE216187 15-subset Spearman similarity heatmap) | `Shared_Inputs/GSE216187_similarity/` | `33_direct_full15_plot.R` |
| `S05ef` | Fig. S5E–F (GSE222633 CD8 UMAP & canonical marker dot plot) | `Shared_Inputs/GSE222633_CD8/` | `21_plot_FCD8.R` |
| `S05f` | Fig. S5F (GSE222633 CD8 revised marker dot plot) | `Shared_Inputs/GSE222633_CD8/` | `25_marker_dotplot_v2.R` |
| `S05g` | Fig. S5G (GSE222633 longitudinal Teff expression & functional scores) | `Shared_Inputs/GSE222633_longitudinal/` | `06_redraw_all_figures.R` |

## 3. Primary Server Extraction Scripts

- [01_export_Tcell_plot_inputs.R](01_export_Tcell_plot_inputs.R) and [02_export_PBMC_plot_inputs.R](02_export_PBMC_plot_inputs.R) were executed directly against the verified server objects (`scRNA224.rds`: 175,531 cells / 23 samples; `Tcell.rds`: 72,884 cells / 15 subsets; `CD8cds.rds`; `CD4cds.rds`). No clustering, annotation, statistical test, or trajectory model was refitted.
- Full multi-gigabyte object locations and public GEO accessions are listed in [00_External_Data_Dependencies.tsv](00_External_Data_Dependencies.tsv).

## 4. Transparent Disclosure of Archived Source Boundaries

1. **Interactive original scripts & `.source.txt` fragments**: Archived R/Rmd scripts retain original server paths and interactive session variable names for provenance. Four historical files containing incomplete or mixed-language notes are explicitly stored as `.source.txt` (`01_gene_set_scores.source.txt`, `02_functional_scores_and_composition.source.txt`, `01_monocle_analysis.source.txt`, `01_CDR3_motifs.source.txt`) rather than standalone runnable programs.
2. **Unlocated external/design sources**:
   - `Fig. 1A`: Editable design source for the experimental workflow schematic was not archived on server 252.
   - `Fig. 4D`: Generated via the web-based Kaplan-Meier Plotter platform (TCGA + GSE12417, GSE37642, GSE6891); exact patient-level exported survival tables were not stored locally.
   - `Fig. S3E–F`: The exact `NeoTCR` scoring script and precomputed score column were not present in the archived `Tcell.rds` metadata or inspected scripts.
   - `GSE222633 FCD8.rds`: Referenced by historical annotation code but absent at its stated top-level server path; all derived `FCD8_*` annotation and plotting tables are included in `Shared_Inputs/GSE222633_CD8/`, and the available `CD8.rds` is documented separately without assuming identity.
3. **Submitted composite figures**: All 11 `Combined/` directories retain the exact original EMF graphics embedded in the submitted manuscript together with verified single-page PDF and 300-DPI PNG renderings.
