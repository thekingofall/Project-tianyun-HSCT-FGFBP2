library(data.table)
library(ComplexHeatmap)
library(circlize)
library(grid)
path <- '/data5/tiany/GSE216187/GSE216187_effector_DE_heatmap_v20261029_013019'
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
d <- fread(file.path(path,'data/plot_ready/heatmap_row_zscores_markers16_v4.tsv'))
z <- as.matrix(d[,-1,with=FALSE]);rownames(z)<-d$gene
patients <- fread(file.path(path,'data/plot_ready/patient_annotation.tsv'))
patients <- patients[match(colnames(z),orig)]
genes_CR <- c('KLRG1', 'HLA-DRB5', 'IKZF1', 'GZMH', 'CX3CR1', 'ZNF683', 'FGFBP2', 'KLRD1', 'TRGC2')
genes_REL <- c('CXCR4', 'DUSP2', 'NR4A2', 'TNFAIP3', 'FTH1', 'IFI6', 'CD27')
left_mark <- rowAnnotation(
  Genes=anno_mark(at=match(genes_REL,rownames(z)),labels=genes_REL,side='left',
    labels_gp=gpar(fontfamily='Arial',fontsize=13),link_gp=gpar(col='black',lwd=1),
    padding=unit(.15,'mm')),show_annotation_name=FALSE)
right_mark <- rowAnnotation(
  Genes=anno_mark(at=match(genes_CR,rownames(z)),labels=genes_CR,side='right',
    labels_gp=gpar(fontfamily='Arial',fontsize=13),link_gp=gpar(col='black',lwd=1),
    padding=unit(.15,'mm')),show_annotation_name=FALSE)
top <- HeatmapAnnotation(
  Group=factor(patients$group,levels=c('REL','CR')),
  col=list(Group=c(REL='#61749B',CR='#B15A50')),
  show_annotation_name=FALSE,simple_anno_size=unit(2,'mm'),
  annotation_legend_param=list(Group=list(title='Outcome',direction='horizontal',nrow=1,
    title_gp=gpar(fontfamily='Arial',fontsize=12),
    labels_gp=gpar(fontfamily='Arial',fontsize=11))))

# Orient the patient dendrogram with REL toward the left, retaining its distances.
patient_tree <- as.dendrogram(hclust(dist(t(z)),method='complete'))
patient_tree <- reorder(patient_tree,wts=ifelse(patients$group=='REL',0,1),agglo.FUN=mean)
ht <- Heatmap(z,name='Z-score',col=colorRamp2(c(-2,0,2),c('#5C53B6','#F7F7F7','#F16643')),
  cluster_rows=TRUE,clustering_distance_rows='euclidean',clustering_method_rows='complete',
  cluster_columns=patient_tree,column_dend_reorder=FALSE,
  show_row_names=FALSE,show_row_dend=FALSE,
  column_names_side='top',column_names_rot=90,
  column_names_gp=gpar(fontfamily='Arial',fontsize=13),
  column_dend_height=unit(8,'mm'),column_dend_gp=gpar(col='black',lwd=1.2),
  top_annotation=top,left_annotation=left_mark,right_annotation=right_mark,
  rect_gp=gpar(col=NA),border=FALSE,
  width=unit(65,'mm'),height=unit(45,'mm'),use_raster=FALSE,
  heatmap_legend_param=list(at=c(-2,-1,0,1,2),legend_width=unit(45,'mm'),direction='horizontal',
    title_gp=gpar(fontfamily='Arial',fontsize=12),labels_gp=gpar(fontfamily='Arial',fontsize=11)))

dir.create(file.path(path,'figures/png'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(path,'figures/pdf'),recursive=TRUE,showWarnings=FALSE)
png_file<-file.path(path,'figures/png/Effector_CR_REL_patient_line_heatmap_markers16_compact_v4.png')
pdf_file<-file.path(path,'figures/pdf/Effector_CR_REL_patient_line_heatmap_markers16_compact_v4.pdf')
stopifnot(!file.exists(png_file),!file.exists(pdf_file))
grDevices::png(png_file,width=7.6,height=3.7666666667,units='in',res=300,type='cairo',bg='white')
drawn <- draw(ht,heatmap_legend_side='bottom',annotation_legend_side='bottom',
              merge_legends=TRUE,padding=unit(c(3,4,3,4),'mm'))
dev.off()
grDevices::cairo_pdf(pdf_file,width=7.6,height=3.7666666667,family='Arial',bg='white')
draw(ht,heatmap_legend_side='bottom',annotation_legend_side='bottom',
     merge_legends=TRUE,padding=unit(c(3,4,3,4),'mm'))
dev.off()
row_id<-row_order(drawn);col_id<-column_order(drawn)
fwrite(data.table(gene=rownames(z)[row_id],z[row_id,col_id,drop=FALSE]),
       file.path(path,'data/plot_ready/heatmap_display_order_markers16_compact_v4.tsv'),sep='\t')
cat('MARKERS16_COMPACT_HEATMAP_READY\n')
