# HSCT T cell remodeling and FGFBP2

This study profiles peripheral immune reconstitution after allogeneic hematopoietic stem cell transplantation using longitudinal single-cell RNA and TCR sequencing. The manuscript describes an FGFBP2-marked effector T cell state and examines its association with clinical response in independent AML cohorts.

The repository follows the supplied November 12 manuscript: **5 main figures, 6 supplementary figures and 71 panels**. The primary PBMC object is `scRNA224.rds` (175,531 cells from 23 samples); the saved T cell object contains 72,884 cells and 15 annotated subsets. The combined manuscript figures remain the reference for the submitted presentation.

[Manuscript](docs/manuscript/Manuscript-11.12.docx) · [Supplementary figures](docs/manuscript/Supplementary_Figures-11.12.docx) · [Reviewer responses](docs/manuscript/Responses-11.12.pdf) · [Panel index](00_Settings/00_Panel_Code_Index.tsv) · [Data and execution](00_Global_Data/README.md)

### Fig. 1 | Longitudinal single-cell atlas of peripheral blood immune reconstitution post-allo-HSCT.

![M01](Main_Figure_01/Combined/M01.png)

[PDF](Main_Figure_01/Combined/M01.pdf) · [Panels, code and data](Main_Figure_01/README.md)

Early reconstitution shifts toward myeloid cells, followed by numerical T cell recovery. The post-transplant T cell distribution remains distinct from healthy donors.

### Fig. 2 | Remodelling of T cell functional states and the TCR repertoire post-allo-HSCT.

![M02](Main_Figure_02/Combined/M02.png)

[PDF](Main_Figure_02/Combined/M02.pdf) · [Panels, code and data](Main_Figure_02/README.md)

T cell recovery is biased toward effector and exhausted states. Functional activation increases while the TCR repertoire remains restricted.

### Fig. 3 | Trajectory analysis identifies distinct state transitions in CD8 and CD4 T cells.

![M03](Main_Figure_03/Combined/M03.png)

[PDF](Main_Figure_03/Combined/M03.pdf) · [Panels, code and data](Main_Figure_03/README.md)

Trajectory and clonal analyses distinguish effector-rich states from naive and exhausted compartments in CD8 and CD4 T cells.

### Fig. 4 | Transcriptional programmes linked to effector-state remodelling and clinical response.

![M04](Main_Figure_04/Combined/M04.png)

[PDF](Main_Figure_04/Combined/M04.pdf) · [Panels, code and data](Main_Figure_04/README.md)

Effector-state transcriptional programs identify FGFBP2 as a candidate marker. External cohorts link its expression to favorable clinical response; these observations support association rather than a causal treatment effect.

### Fig. 5 | Clonal architecture and T cell receptor convergence across CD8 effector states.

![M05](Main_Figure_05/Combined/M05.png)

[PDF](Main_Figure_05/Combined/M05.pdf) · [Panels, code and data](Main_Figure_05/README.md)

Clonal overlap, receptor usage and CDR3 motifs characterize continuity and differences between CD8 effector and exhausted states.

### Supplementary Figure 1 | Overview of sample distribution and immune lineage composition following HSCT.

![S01](Supplementary_Figures/Supplementary_Figure_01/Combined/S01.png)

[PDF](Supplementary_Figures/Supplementary_Figure_01/Combined/S01.pdf) · [Panels, code and data](Supplementary_Figures/Supplementary_Figure_01/README.md)

Sample contributions and lineage proportions provide the compositional context for the longitudinal atlas.

### Supplementary Figure 2 | Dynamic changes in T cell composition and spatial distribution post-allo-HSCT.

![S02](Supplementary_Figures/Supplementary_Figure_02/Combined/S02.png)

[PDF](Supplementary_Figures/Supplementary_Figure_02/Combined/S02.pdf) · [Panels, code and data](Supplementary_Figures/Supplementary_Figure_02/README.md)

The CD4/CD8 balance, subset composition and functional scores show persistent changes during reconstitution.

### Supplementary Figure 3 | Compositional analysis of T cell developmental states and external validation of transcriptomic signatures.

![S03](Supplementary_Figures/Supplementary_Figure_03/Combined/S03.png)

[PDF](Supplementary_Figures/Supplementary_Figure_03/Combined/S03.pdf) · [Panels, code and data](Supplementary_Figures/Supplementary_Figure_03/README.md)

State-specific composition and NeoTCR scores further characterize the inferred T cell states.

### Supplementary Figure 4 | T cell landscape and enrichment analysis in the independent clinical validation cohort.

![S04](Supplementary_Figures/Supplementary_Figure_04/Combined/S04.png)

[PDF](Supplementary_Figures/Supplementary_Figure_04/Combined/S04.pdf) · [Panels, code and data](Supplementary_Figures/Supplementary_Figure_04/README.md)

The GSE223844 cohort provides an independent view of treatment response, T cell composition and FGFBP2 expression.

### Supplementary Figure 5. Independent validation of the FGFBP2-associated effector T cell phenotype in AML allo-HSCT cohorts.

![S05](Supplementary_Figures/Supplementary_Figure_05/Combined/S05.png)

[PDF](Supplementary_Figures/Supplementary_Figure_05/Combined/S05.pdf) · [Panels, code and data](Supplementary_Figures/Supplementary_Figure_05/README.md)

GSE216187 and GSE222633 extend the effector-state comparison to independent transplant cohorts.

### Supplementary Figure 6 | Clonal lineage connectivity and expression profiles of immune checkpoints across T cell states.

![S06](Supplementary_Figures/Supplementary_Figure_06/Combined/S06.png)

[PDF](Supplementary_Figures/Supplementary_Figure_06/Combined/S06.pdf) · [Panels, code and data](Supplementary_Figures/Supplementary_Figure_06/README.md)

Receptor gene usage and inhibitory-marker expression complement the clonal and transcriptional comparisons.

## Data and execution

Cell annotations, UMAP and trajectory coordinates, functional scores, selected normalized gene expression and revision plotting tables are included. Multi-gigabyte Seurat objects and raw sequencing data remain on server 252; [external dependencies](00_Settings/00_External_Data_Dependencies.tsv) records the exact available object paths and cohort accessions. [File provenance](00_Settings/01_File_Provenance.tsv) identifies copied source files.

From the repository root, prepare a saved heatmap in a new folder:

```bash
python 00_Settings/run_saved_plot.py M04h --output /path/to/new-output --prepare-only
```

Omit `--prepare-only` to execute with Rscript. The other supported table entries are `S05ef` and `S05f`. R requires the packages imported by the chosen script and Arial for the original typography. Input preparation was tested; these plotting entries have not been rerun during this organization task.

To regenerate primary small inputs on server 252 from the existing fitted objects:

```bash
Rscript 00_Settings/01_export_Tcell_plot_inputs.R /data5/tiany/HSCT/Result /path/to/new-input-folder
Rscript 00_Settings/02_export_PBMC_plot_inputs.R /data5/tiany/HSCT/Result /path/to/new-input-folder
```

The extraction was executed against the server objects without refitting or changing annotations. [Execution notes](00_Settings/README.md) distinguish saved-table entries, interactive original code and source fragments. The diagram design file, exact Kaplan-Meier Plotter patient-level input and some standalone final panel exports were not located. Their submitted images are retained; no replacement data were constructed.
