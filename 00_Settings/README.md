# Inputs and execution

The panel index links all 71 submitted panels to available code, data and source notes. Original interactive code retains its source paths and object names for provenance; it should be read section by section with the matching loaded objects, not sourced as an end-to-end pipeline. Source paths and full objects are listed in the external dependency table.

Four historical files contain incomplete or mixed-language fragments and are stored as `.source.txt`: gene-set scoring notes, functional/composition notes, Monocle notes, and R/Python CDR3 motif notes. They are not advertised as standalone R programs. The remaining original R files and extracted R Markdown chunks passed R parsing; parsing does not establish end-to-end execution.

`run_saved_plot.py` resolves the repository from its own location, copies small input tables to a new output folder, substitutes runtime paths and runs the archived table script there. It supports M04h (heatmap), S05ef (CD8 UMAP and marker dot plot), and S05f (revised marker dot plot). Its input preparation was verified from multiple working directories; plotting was not rerun. Other archived figures include historical plot types and are retained as supplied.

The two export R scripts were run against the primary server objects. Their outputs preserve cell identities and stored values. No clustering, statistical testing or trajectory fitting was rerun. Software details from this extraction are recorded in R_sessionInfo.txt.

Known source gaps: the editable M01a schematic; exact patient-level input used by the Kaplan-Meier Plotter result in M04d; and standalone exports for several primary panels. Original interactive code contains reused working variables, so an exact run command is supplied only where a saved-table entry is supported. Final composite figures come directly from the submitted Word documents.

The historical GSE222633 annotation script names FCD8.rds as an output, but that file was not present at its stated top-level server path during this inventory. Its saved annotation and plotting tables are included; the available CD8.rds is listed separately without assuming it is identical.

The exact NeoTCR scoring implementation and score column for S03e/f were not located in the saved Tcell object or inspected scripts; the submitted panels are preserved.
