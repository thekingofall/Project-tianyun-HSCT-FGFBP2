# Supplementary Figure 5. Independent validation of the FGFBP2-associated effector T cell phenotype in AML allo-HSCT cohorts.

GSE216187 and GSE222633 extend the effector-state comparison to independent transplant cohorts.

![S05](Combined/S05.png)

[PDF](Combined/S05.pdf) · [Original EMF](Combined/S05.emf)

(A) UMAP visualization of annotated bone marrow cell populations from six patients with AML approximately 100 days after allo-HSCT, including three patients in CR and three patients who subsequently experienced early relapse (REL).

(B) Heatmap showing the expression of representative marker genes across the annotated CD8⁺ T cell, CD4⁺ T cell, and hematopoietic stem and progenitor cell (HSPC) populations.

(C) Spearman correlation analysis assessing the transcriptional similarity between CD8⁺ T cell subsets in the external cohort and the S2_CD8_Teff population identified in our study.

(D) UMAP feature plots showing the expression distributions of ADGRG1 and FGFBP2 in the CR and REL groups. Panels A–D were generated from the external validation cohort (GEO: GSE216187).

(E) UMAP visualization of annotated CD8⁺ T cell subsets, including naïve, Tcm, Tem, Teff, and cycling T cells, as well as Trm-like cells.

(F) Dot plot showing the expression of representative marker genes across the annotated CD8⁺ T cell subsets. Dot size indicates the percentage of cells expressing each gene, and color represents the average expression level.

(G) Violin plots comparing FGFBP2 expression in CD8⁺ Teff cells between the CR and REL groups across the indicated sampling time points (post-G-CSF mobilization, and days +90 and +180 post-allo-HSCT). Red and blue indicate the CR and REL groups, respectively. Panels E–G were generated from the external validation cohort (GEO: GSE222633).

## Panel code and data

[A](a_S05a/README.md) · [B](b_S05b/README.md) · [C](c_S05c/README.md) · [D](d_S05d/README.md) · [E](e_S05e/README.md) · [F](f_S05f/README.md) · [G](g_S05g/README.md)

## Figure source and execution

Combined contains the exact EMF embedded in the supplied manuscript and raster PNG/PDF previews. Server panel exports are retained where located. Panel pages identify shared inputs and original scripts. Historical interactive scripts require the documented session objects; source fragments ending in `.source.txt` are not standalone programs.

Table-based CD8 UMAP and marker entries from the repository root: `python 00_Settings/run_saved_plot.py S05ef --output /path/to/new-output`.
