#!/usr/bin/env Rscript
# Portable saved-table reproduction for primary cohort panels (Fig. 1-3, S1-S3, S6E)
# Reads ONLY the compact exported tables under 00_Global_Data/Shared_Inputs/Primary/
# and 00_Settings/global_palette.tsv; no multi-gigabyte Seurat RDS files required.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript 03_reproduce_primary_panels.R <module: M01|M02|M03|S06e|all> <repo_root> <output_dir>")
}
module <- args[1]
repo_root <- normalizePath(args[2], mustWork = TRUE)
out_dir <- args[3]

primary_dir <- file.path(repo_root, "00_Global_Data", "Shared_Inputs", "Primary")
palette_file <- file.path(repo_root, "00_Settings", "global_palette.tsv")
pal_dt <- fread(palette_file)

get_pal <- function(pid) {
  sub <- pal_dt[palette_id == pid][order(order)]
  setNames(sub$hex, sub$label)
}

dir.create(file.path(out_dir, "figures", "png"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures", "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

base_family <- if (Sys.info()[["sysname"]] == "Linux") "sans" else "Arial"
theme_paper <- function(base_size = 13) {
  theme_classic(base_family = base_family, base_size = base_size) +
    theme(
      axis.line = element_line(colour = "black", linewidth = 0.7),
      axis.ticks = element_line(colour = "black", linewidth = 0.6),
      axis.text = element_text(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.key.size = unit(4.5, "mm"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

manifest <- list()
save_pair <- function(p, id, w, h, src) {
  png_path <- file.path(out_dir, "figures", "png", paste0(id, ".png"))
  pdf_path <- file.path(out_dir, "figures", "pdf", paste0(id, ".pdf"))
  png_type <- if (isTRUE(capabilities("cairo"))) "cairo" else "Xlib"
  grDevices::png(png_path, width = w, height = h, units = "in", res = 300, type = png_type, bg = "white")
  print(p)
  grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = w, height = h, bg = "white")
  print(p)
  grDevices::dev.off()
  manifest[[id]] <<- data.table(panel_id = id, png = png_path, pdf = pdf_path, width_in = w, height_in = h, source_input = src)
}

run_M01 <- function() {
  meta <- fread(file.path(primary_dir, "PBMC_metadata.csv.gz"))
  umap <- fread(file.path(primary_dir, "PBMC_UMAP.csv.gz"))
  expr <- fread(file.path(primary_dir, "PBMC_marker_logexpression.csv.gz"))
  df <- merge(meta[, .(cell_id, orig.ident, Celltype, condition)], umap, by = "cell_id", sort = FALSE)
  pbmc_cols <- get_pal("Primary_PBMC")
  stage_order <- c("HD", "Pre-HSCT", "Post-1M", "Post-3M", "Post-6M")
  df[, Celltype := factor(Celltype, levels = names(pbmc_cols))]
  df[, condition := factor(condition, levels = stage_order)]

  # M01b: PBMC UMAP
  set.seed(2026)
  df_plot <- df[sample(.N)]
  p_m01b <- ggplot(df_plot, aes(UMAP_1, UMAP_2, colour = Celltype)) +
    geom_point(size = 0.15, alpha = 0.75, stroke = 0) +
    scale_colour_manual(values = pbmc_cols, name = "Lineage") +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "Fig. 1B | Integrated PBMC landscape (n = 175,531 cells)") +
    theme_paper()
  save_pair(p_m01b, "M01b_PBMC_UMAP", 8.2, 6.2, "PBMC_metadata.csv.gz;PBMC_UMAP.csv.gz")

  # M01c: Canonical marker dot plot
  markers <- c("CD34","SPINK2","PROM1","CD3E","CD3D","CD3G","GNLY","NKG7","CX3CR1",
               "CD79A","CD19","MS4A1","CD14","VCAN","LYZ","FCGR3B","CMTM2","S100A8",
               "CD1C","FCER1A","ENHO","IL3RA","CLEC4C","JCHAIN","PPBP","PF4","GP9")
  markers <- intersect(markers, names(expr))
  ex_m <- merge(meta[, .(cell_id, Celltype)], expr[, c("cell_id", markers), with = FALSE], by = "cell_id", sort = FALSE)
  long_m <- melt(ex_m, id.vars = c("cell_id", "Celltype"), variable.name = "Gene", value.name = "Expr")
  dot_m <- long_m[, .(avg_exp = mean(Expr), pct_exp = mean(Expr > 0) * 100), by = .(Celltype, Gene)]
  dot_m[, scaled_exp := as.numeric(scale(avg_exp)), by = Gene]
  dot_m[, Celltype := factor(Celltype, levels = rev(names(pbmc_cols)))]
  dot_m[, Gene := factor(Gene, levels = markers)]
  fwrite(dot_m, file.path(out_dir, "tables", "M01c_PBMC_marker_dotplot.tsv"), sep = "\t")
  p_m01c <- ggplot(dot_m, aes(Gene, Celltype, size = pct_exp, colour = scaled_exp)) +
    geom_point() +
    scale_size_area(max_size = 5.5, limits = c(0, 100), name = "Expressed (%)") +
    scale_colour_gradient2(low = "#326A82", mid = "#F7F7F7", high = "#B52635", midpoint = 0, limits = c(-2, 2), oob = squish, name = "Scaled avg.") +
    labs(x = NULL, y = NULL, title = "Fig. 1C | Canonical PBMC lineage markers") +
    theme_paper(11) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1), legend.position = "bottom")
  save_pair(p_m01c, "M01c_PBMC_Marker_DotPlot", 10.5, 4.8, "PBMC_metadata.csv.gz;PBMC_marker_logexpression.csv.gz")

  # M01d & S01b: Longitudinal lineage composition (stacked area & stacked bar)
  comp <- df[, .N, by = .(condition, Celltype)]
  comp[, prop := N / sum(N) * 100, by = condition]
  comp[, stage_x := as.integer(condition)]
  fwrite(comp, file.path(out_dir, "tables", "M01d_PBMC_lineage_composition.tsv"), sep = "\t")
  p_m01d <- ggplot(comp, aes(stage_x, prop, fill = Celltype)) +
    geom_area(colour = "white", linewidth = 0.4) +
    scale_fill_manual(values = pbmc_cols, name = "Lineage") +
    scale_x_continuous(breaks = seq_along(stage_order), labels = stage_order, expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion (%)", title = "Fig. 1D | Longitudinal immune lineage composition") +
    theme_paper()
  save_pair(p_m01d, "M01d_PBMC_Composition_Area", 7.8, 5.0, "PBMC_metadata.csv.gz")

  # M01e: Timepoint UMAP density
  p_m01e <- ggplot(df_plot, aes(UMAP_1, UMAP_2)) +
    stat_density_2d(aes(fill = after_stat(ndensity)), geom = "raster", contour = FALSE) +
    scale_fill_distiller(palette = "Blues", direction = 1, name = "Relative\ndensity") +
    facet_wrap(~condition, nrow = 1) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "Fig. 1E | PBMC density redistribution across sampling stages") +
    theme_paper(11) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))
  save_pair(p_m01e, "M01e_PBMC_Stage_Density", 12.5, 3.2, "PBMC_metadata.csv.gz;PBMC_UMAP.csv.gz")
}

run_M02 <- function() {
  meta <- fread(file.path(primary_dir, "Tcell_metadata.csv.gz"))
  expr <- fread(file.path(primary_dir, "Tcell_marker_logexpression.csv.gz"))
  t_cols <- get_pal("Primary_Tcell")
  tp_cols <- get_pal("Primary_Timepoint")
  names(tp_cols)[names(tp_cols) == "Pre-HSCT"] <- "Pre"
  stage_order <- c("Pre", "HD", "Post-1M", "Post-3M", "Post-6M")

  meta[, Celltype := factor(Celltype, levels = names(t_cols))]
  meta[, condition := factor(condition, levels = stage_order)]

  # M02a: T cell UMAP
  set.seed(2026)
  m_plot <- meta[sample(.N)]
  p_m02a <- ggplot(m_plot, aes(UMAP_1, UMAP_2, colour = Celltype)) +
    geom_point(size = 0.18, alpha = 0.8, stroke = 0) +
    scale_colour_manual(values = t_cols, name = "T cell subset") +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2", title = "Fig. 2A | T cell sub-populations (n = 72,884 cells)") +
    theme_paper()
  save_pair(p_m02a, "M02a_Tcell_UMAP", 8.8, 6.4, "Tcell_metadata.csv.gz")

  # M02b: T cell canonical markers
  t_markers <- c("CD3D","CD4","CD8A","CD8B","IL7R","CCR7","TCF7","SELL","FOXP3","IL2RA",
                 "GZMK","GNLY","PRF1","GZMH","GZMB","FGFBP2","ADGRG1","LAG3","TIGIT","PDCD1",
                 "ZNF683","KLRB1","SLC4A10","TRDC","TRGC1","TRGC2","MKI67","TOP2A")
  t_markers <- intersect(t_markers, names(expr))
  ex_t <- merge(meta[, .(cell_id, Celltype)], expr[, c("cell_id", t_markers), with = FALSE], by = "cell_id", sort = FALSE)
  long_t <- melt(ex_t, id.vars = c("cell_id", "Celltype"), variable.name = "Gene", value.name = "Expr")
  dot_t <- long_t[, .(avg_exp = mean(Expr), pct_exp = mean(Expr > 0) * 100), by = .(Celltype, Gene)]
  dot_t[, scaled_exp := as.numeric(scale(avg_exp)), by = Gene]
  dot_t[, Celltype := factor(Celltype, levels = rev(names(t_cols)))]
  dot_t[, Gene := factor(Gene, levels = t_markers)]
  fwrite(dot_t, file.path(out_dir, "tables", "M02b_Tcell_marker_dotplot.tsv"), sep = "\t")
  p_m02b <- ggplot(dot_t, aes(Gene, Celltype, size = pct_exp, colour = scaled_exp)) +
    geom_point() +
    scale_size_area(max_size = 5, limits = c(0, 100), name = "Expressed (%)") +
    scale_colour_gradient2(low = "#326A82", mid = "#F7F7F7", high = "#B52635", midpoint = 0, limits = c(-2, 2), oob = squish, name = "Scaled avg.") +
    labs(x = NULL, y = NULL, title = "Fig. 2B | Lineage and functional marker expression across T cell subsets") +
    theme_paper(11) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1), legend.position = "bottom")
  save_pair(p_m02b, "M02b_Tcell_Marker_DotPlot", 11.2, 6.2, "Tcell_metadata.csv.gz;Tcell_marker_logexpression.csv.gz")

  # M02c: Temporal T cell subset composition
  stage_m02c <- c("HD", "Pre", "Post-1M", "Post-3M", "Post-6M")
  comp <- meta[, .N, by = .(condition, Celltype)]
  comp[, condition := factor(condition, levels = stage_m02c)]
  comp[, prop := N / sum(N) * 100, by = condition]
  comp[, stage_x := as.integer(condition)]
  fwrite(comp, file.path(out_dir, "tables", "M02c_Tcell_subset_composition.tsv"), sep = "\t")
  p_m02c <- ggplot(comp, aes(stage_x, prop, fill = Celltype)) +
    geom_area(colour = "white", linewidth = 0.35) +
    scale_fill_manual(values = t_cols, name = "T cell subset") +
    scale_x_continuous(breaks = seq_along(stage_m02c), labels = stage_m02c, expand = c(0.01, 0.01)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion (%)", title = "Fig. 2C | Temporal dynamics of T cell subset composition") +
    theme_paper()
  save_pair(p_m02c, "M02c_Tcell_Composition_Area", 8.5, 5.4, "Tcell_metadata.csv.gz")

  # M02d & S02d: Functional scores across stages
  score_cols <- intersect(c("activation1", "cytotoxicity1", "texterm1", "TCRsign1", "cycle1"), names(meta))
  score_labels <- c(activation1 = "Activation", cytotoxicity1 = "Cytotoxicity",
                    texterm1 = "Terminal exhaustion", TCRsign1 = "TCR signaling", cycle1 = "Proliferation")
  sc_long <- melt(meta[, c("condition", score_cols), with = FALSE], id.vars = "condition",
                  variable.name = "Program", value.name = "Score")
  sc_long[, Program := factor(score_labels[as.character(Program)], levels = unname(score_labels[score_cols]))]
  sc_summary <- sc_long[, .(n = .N, median = median(Score, na.rm = TRUE),
                            q25 = quantile(Score, 0.25, na.rm = TRUE),
                            q75 = quantile(Score, 0.75, na.rm = TRUE)), by = .(Program, condition)]
  fwrite(sc_summary, file.path(out_dir, "tables", "M02d_S02d_functional_score_summary.tsv"), sep = "\t")
  p_m02d <- ggplot(sc_long, aes(condition, Score, fill = condition)) +
    geom_violin(scale = "width", colour = "black", linewidth = 0.4, trim = TRUE) +
    geom_boxplot(width = 0.14, outlier.shape = NA, colour = "black", fill = "white", linewidth = 0.35) +
    facet_wrap(~Program, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = tp_cols, name = "Stage") +
    labs(x = NULL, y = "Module score", title = "Fig. 2D & Fig. S2D | Functional scores across transplant stages") +
    theme_paper(11) +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "none")
  save_pair(p_m02d, "M02d_S02d_Functional_Scores", 13.5, 4.2, "Tcell_metadata.csv.gz")
}

run_M03 <- function() {
  t_cols <- get_pal("Primary_Tcell")
  st_cols <- get_pal("Primary_Trajectory_State")
  stage_cols <- c(HD = "#0070C0", Pre = "#7030A0", Post = "#E25A4C")

  for (lin in c("CD8", "CD4")) {
    md <- fread(file.path(primary_dir, paste0(lin, "cds_metadata.csv.gz")))
    co <- fread(file.path(primary_dir, paste0(lin, "cds_trajectory_coordinates.csv.gz")))
    mp <- fread(file.path(primary_dir, paste0(lin, "_raw_state_to_saved_annotation.csv")))
    df <- merge(md[, .(cell_id, Celltype, condition, Pseudotime, State)], co, by = "cell_id", sort = FALSE)
    # Recover S1/S2/S3 from exact raw-to-saved mapping table
    map_vec <- setNames(as.character(mp$Saved_Tcell_state), as.character(mp$Raw_Monocle_state))
    df[, SavedState := factor(map_vec[as.character(State)], levels = c("S1", "S2", "S3"))]
    df[, StageGroup := factor(ifelse(condition %in% c("Post-1M", "Post-3M", "Post-6M"), "Post", condition),
                              levels = c("HD", "Pre", "Post"))]
    sub_cols <- t_cols[intersect(names(t_cols), unique(df$Celltype))]
    df[, Celltype := factor(Celltype, levels = names(sub_cols))]

    prefix <- if (lin == "CD8") "M03abc" else "M03ghi"
    p_ct <- ggplot(df, aes(V1, V2, colour = Celltype)) +
      geom_point(size = 0.25, alpha = 0.8, stroke = 0) +
      scale_colour_manual(values = sub_cols, name = paste(lin, "subset")) +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(x = "Component 1", y = "Component 2", title = paste0(lin, " trajectory: Cell type (n = ", nrow(df), ")")) +
      theme_paper(11)
    save_pair(p_ct, paste0(prefix, "_", lin, "_Trajectory_Celltype"), 7.2, 5.2, paste0(lin, "cds_metadata.csv.gz;", lin, "cds_trajectory_coordinates.csv.gz"))

    p_st <- ggplot(df[!is.na(SavedState)], aes(V1, V2, colour = SavedState)) +
      geom_point(size = 0.25, alpha = 0.8, stroke = 0) +
      scale_colour_manual(values = st_cols, name = "Inferred state") +
      guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
      labs(x = "Component 1", y = "Component 2", title = paste0(lin, " trajectory: Inferred states (S1-S3)")) +
      theme_paper(11)
    save_pair(p_st, paste0(prefix, "_", lin, "_Trajectory_State"), 6.8, 5.2, paste0(lin, "cds_metadata.csv.gz;", lin, "_raw_state_to_saved_annotation.csv"))

    # State composition by Celltype (M03d / M03j) and by StageGroup (M03e / M03k / S03a-d)
    st_ct <- df[!is.na(SavedState), .N, by = .(SavedState, Celltype)]
    st_ct[, prop := N / sum(N) * 100, by = SavedState]
    fwrite(st_ct, file.path(out_dir, "tables", paste0("M03_", lin, "_state_subset_composition.tsv")), sep = "\t")
    p_comp <- ggplot(st_ct, aes(SavedState, prop, fill = Celltype)) +
      geom_col(width = 0.68, colour = "white", linewidth = 0.3) +
      scale_fill_manual(values = sub_cols, name = paste(lin, "subset")) +
      labs(x = "Inferred state", y = "Proportion (%)", title = paste0(lin, " subset distribution across states")) +
      theme_paper(11)
    save_pair(p_comp, paste0("M03_", lin, "_State_Composition"), 6.2, 4.8, paste0(lin, "cds_metadata.csv.gz"))

    st_cond <- df[!is.na(SavedState), .N, by = .(SavedState, Celltype, StageGroup)]
    st_cond[, prop := N / sum(N) * 100, by = .(SavedState, Celltype)]
    fwrite(st_cond, file.path(out_dir, "tables", paste0("M03_S03_", lin, "_state_stage_composition.tsv")), sep = "\t")
    p_cond <- ggplot(st_cond, aes(Celltype, prop, fill = StageGroup)) +
      geom_col(width = 0.7, colour = "white", linewidth = 0.25) +
      facet_wrap(~SavedState, nrow = 1, scales = "free_x") +
      scale_fill_manual(values = stage_cols, name = "Stage") +
      labs(x = NULL, y = "Proportion (%)", title = paste0(lin, " stage contribution across S1-S3 subsets")) +
      theme_paper(10.5) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_pair(p_cond, paste0("M03_S03_", lin, "_State_Stage_Breakdown"), 10.5, 4.6, paste0(lin, "cds_metadata.csv.gz"))
  }
}

run_S06e <- function() {
  meta <- fread(file.path(primary_dir, "Tcell_metadata.csv.gz"))
  expr <- fread(file.path(primary_dir, "Tcell_marker_logexpression.csv.gz"))
  genes <- intersect(c("LAG3", "TIGIT", "PDCD1", "CTLA4", "HAVCR2", "TOX"), names(expr))
  target_states <- c("S1_CD8_Naive", "S2_CD8_Teff", "S3_CD8_Teff", "S3_CD8_Tex",
                     "S1_CD4_Naive", "S2_CD4_Teff", "S3_CD4_Treg")
  sub <- merge(meta[Celltype2 %in% target_states, .(cell_id, Celltype2)],
               expr[, c("cell_id", genes), with = FALSE], by = "cell_id", sort = FALSE)
  long_g <- melt(sub, id.vars = c("cell_id", "Celltype2"), variable.name = "Gene", value.name = "Expr")
  pos <- long_g[, .(n_cells = .N, pct_pos = mean(Expr > 0) * 100, mean_expr = mean(Expr)), by = .(Celltype2, Gene)]
  pos[, Celltype2 := factor(Celltype2, levels = target_states)]
  pos[, Gene := factor(Gene, levels = genes)]
  fwrite(pos, file.path(out_dir, "tables", "S06e_inhibitory_marker_positivity.tsv"), sep = "\t")
  p_s06e <- ggplot(pos, aes(Celltype2, pct_pos, fill = Celltype2)) +
    geom_col(width = 0.72, colour = "black", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.1f%%", pct_pos)), vjust = -0.35, size = 2.8) +
    facet_wrap(~Gene, scales = "free_y", ncol = 3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
    labs(x = NULL, y = "Marker-positive cells (%)", title = "Fig. S6E | Inhibitory receptor and TOX positivity across T cell states") +
    theme_paper(10.5) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
  save_pair(p_s06e, "S06e_Checkpoint_Positivity", 11.2, 6.4, "Tcell_metadata.csv.gz;Tcell_marker_logexpression.csv.gz")
}

targets <- if (module == "all") c("M01", "M02", "M03", "S06e") else module
for (t in targets) {
  if (t == "M01") run_M01()
  else if (t == "M02") run_M02()
  else if (t == "M03") run_M03()
  else if (t == "S06e") run_S06e()
  else stop("Unsupported primary module: ", t)
}

fwrite(rbindlist(manifest), file.path(out_dir, "primary_figure_manifest.tsv"), sep = "\t")
cat("PRIMARY_SAVED_TABLE_PLOTS_READY:", length(manifest), "figures\n")
