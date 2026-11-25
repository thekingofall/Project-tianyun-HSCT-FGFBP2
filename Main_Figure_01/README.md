# Fig. 1 | Longitudinal single-cell atlas of peripheral blood immune reconstitution post-allo-HSCT

> [Back to Repository Overview (EN)](../README.md) · [返回仓库首页 (中文)](../README.zh-CN.md)

![Fig. 1](Combined/M01.png)

[Download PDF](Combined/M01.pdf) · [Original EMF](Combined/M01.emf) · [Primary Shared Inputs](../00_Global_Data/Shared_Inputs/Primary/README.md)

## Key Findings

Longitudinal profiling of 23 peripheral blood mononuclear cell (PBMC) samples from five patients with acute myeloid leukemia (AML) and their paired healthy donors (HD) captures **175,531 high-quality single cells** spanning nine major hematopoietic and immune lineages (`Fig. 1A–C`). Early post-transplant reconstitution (`Post-1M`) is characterized by transient myeloid predominance and a sharp contraction of the T cell compartment (`Fig. 1D–E`). Although total T cell abundance progressively recovers between months 3 and 6 (`Post-3M` to `Post-6M`) to approach HD levels, the transcriptional density centroid of post-transplant T cells remains clearly displaced from both HD and Pre-HSCT baselines, indicating qualitative functional remodeling beyond numerical recovery.

## Panel Descriptions

- **(A)** Experimental design and longitudinal sampling schematic across donor pre-mobilization (`HD`), recipient pre-conditioning (`Pre-HSCT`), and `Post-1M`, `Post-3M`, and `Post-6M` paired scRNA-seq and scTCR-seq.
- **(B)** Integrated UMAP of the global PBMC immune landscape colored by annotated cell lineage (`HSC/MPPs`, `T`, `NK`, `B`, `Monocytes`, `Granulocytes`, `cDC`, `pDC`, and `Platelets`).
- **(C)** Canonical lineage marker dot plot showing scaled average expression and expressing-cell fraction across clusters.
- **(D)** Longitudinal lineage composition across `HD`, `Pre-HSCT`, `Post-1M`, `Post-3M`, and `Post-6M` (stacked area plot).
- **(E)** Stage-stratified UMAP projections and 2D density distributions illustrating dynamic spatial redistribution of PBMCs during reconstitution.

## Panel Directory

[A (`a_M01a`)](a_M01a/README.md) · [B (`b_M01b`)](b_M01b/README.md) · [C (`c_M01c`)](c_M01c/README.md) · [D (`d_M01d`)](d_M01d/README.md) · [E (`e_M01e`)](e_M01e/README.md)

## Code, Inputs & Reproduction

- **Saved-table verification (no large RDS required)**:
  ```bash
  python 00_Settings/run_saved_plot.py M01 --output /path/to/output_M01
  ```
  Generates `M01b` (UMAP), `M01c` (canonical marker dot plot), `M01d` (lineage composition), and `M01e` (stage density maps) directly from `00_Global_Data/Shared_Inputs/Primary/PBMC_*.csv.gz` using locked colors in [GLOBAL_PALETTE.md](../00_Settings/GLOBAL_PALETTE.md).
- **Original archived scripts**: [01_primary_scRNA_TCR_analysis.Rmd](00_Preprocess/01_primary_scRNA_TCR_analysis.Rmd), [01_primary_landscape_and_TCR.R](01_Code/01_primary_landscape_and_TCR.R).
- **Primary server object**: `/data5/tiany/HSCT/Result/scRNA224.rds` (175,531 cells, 23 samples; see [External Data Dependencies](../00_Settings/00_External_Data_Dependencies.tsv)).
- **Known source boundary**: The editable design source for schematic panel `Fig. 1A` was not archived on server 252; the submitted manuscript composite (`Combined/M01.emf`) is retained verbatim.
