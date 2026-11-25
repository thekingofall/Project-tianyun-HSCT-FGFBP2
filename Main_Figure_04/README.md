# Fig. 4 | Transcriptional programmes linked to effector-state remodelling and clinical response.

Effector-state transcriptional programs identify FGFBP2 as a candidate marker. External cohorts link its expression to favorable clinical response; these observations support association rather than a causal treatment effect.

![M04](Combined/M04.png)

[PDF](Combined/M04.pdf) · [Original EMF](Combined/M04.emf)

(A) Heatmaps showing the inferred transcription factor (TF) activity (VIPER) for CD8+ and CD4+ T cells. 

(B) Gene regulatory network (GRN) visualization of differentially active transcriptional regulators (TBX21, EOMES, and PRDM1) and their target differentially expressed genes (DEGs) in CD8+ and CD4+ T cells. 

(C) Volcano plots identifying DEGs between State 2 (S2) and State 3 (S3) Teff cells for both CD8+ and CD4+ lineages. Key effector molecules such as FGFBP2 and ADGRG1 are highlighted. 

(D) Kaplan-Meier survival analysis of an independent AML patient cohort stratified by the expression level of FGFBP2. P values were determined using the log-rank test. 

(E) Three-dimensional summary of the clonal expansion index, tissue distribution index, and state transition index across different CD8+ T cell subsets in patients with complete remission (CR) and non-response (NR).

(F) UMAP projection of CD8⁺ Teff subclusters. Based on their transcriptional profiles, CD8⁺ Teff cells were classified into three subclusters (Teff-1, Teff-2, and Teff-3), and their distributions across the four clinical categories (CR-Pre, CR-Post, NR-Pre, and NR-Post) are shown.

(G) Dot plot illustrating the expression of lineage markers across Teff subclusters. Panels E-G were generated from the external validation cohort (GEO: GSE223844).

(H) Heatmap showing the expression patterns of selected differentially expressed genes in the CD8⁺ T cell compartment between the CR and Relapse (REL) groups. Expression values are displayed as Z-scores (GEO: GSE216187).

(I) Patient-level longitudinal changes in FGFBP2 expression across the same sampling time points. Each line represents an individual donor–recipient pair, with red and blue indicating the CR and REL groups, respectively (GEO: GSE222633).

## Panel code and data

[A](a_M04a/README.md) · [B](b_M04b/README.md) · [C](c_M04c/README.md) · [D](d_M04d/README.md) · [E](e_M04e/README.md) · [F](f_M04f/README.md) · [G](g_M04g/README.md) · [H](h_M04h/README.md) · [I](i_M04i/README.md)

## Figure source and execution

Combined contains the exact EMF embedded in the supplied manuscript and raster PNG/PDF previews. Server panel exports are retained where located. Panel pages identify shared inputs and original scripts. Historical interactive scripts require the documented session objects; source fragments ending in `.source.txt` are not standalone programs.

Table-based heatmap entry from the repository root: `python 00_Settings/run_saved_plot.py M04h --output /path/to/new-output`.
