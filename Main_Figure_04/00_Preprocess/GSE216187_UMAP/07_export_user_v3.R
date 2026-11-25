local({
  oldwd<-getwd();on.exit(setwd(oldwd))
  setwd('/data5/tiany/GSE216187/GSE216187_UMAP_labels_v20261027_204550')
  Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
  Sys.unsetenv('MAO_ARIAL_FONT_DIR')
  png<-'figures/png/GSE216187_UMAP_usercolors_largepoints_v3.png'
  pdf<-'figures/pdf/GSE216187_UMAP_usercolors_largepoints_v3.pdf'
  stopifnot(!file.exists(png),!file.exists(pdf),!file.exists('reports/render_user_v3.json'))
  Bseurat<-get('GSE216187_annotated',envir=.GlobalEnv)
  source('code/06_plain_plot_user_v3.R',local=TRUE)
  stopifnot(nrow(df)==8501L,nrow(label_pos)==24L,!anyNA(df$Celltype))
  ggplot2::ggsave(png,p,width=13.2,height=10.6,units='in',dpi=300,
                 device=grDevices::png,type='cairo',bg='white')
  ggplot2::ggsave(pdf,p,width=13.2,height=10.6,units='in',
                 device=grDevices::cairo_pdf,family='Arial',bg='white')
  data.table::fwrite(label_pos,'data/plot_ready/manual_labels_user_v3.tsv',sep='\t')
  data.table::fwrite(data.frame(Celltype=names(cell_colors),hex=unname(cell_colors)),
                    'data/plot_ready/palette_user_v3.tsv',sep='\t')
  jsonlite::write_json(list(status='complete',cells=nrow(df),labels=nrow(label_pos),
    point_size=1.25,alpha=1,previous_point_size=.72,palette='Exact user-provided HEX values',
    label_placement='Manually reviewed coordinates; leaders anchor to real cells',
    source_object='/data5/tiany/GSE216187/GSE216187_seurat_scRNA_SCT_orig_Celltype_v20261027.rds',
    png=png,pdf=pdf,time=format(Sys.time(),'%FT%T%z')),
    'reports/render_user_v3.json',auto_unbox=TRUE,pretty=TRUE)
  cat('USER_STYLE_UMAP_V3_COMPLETE\n')
})
