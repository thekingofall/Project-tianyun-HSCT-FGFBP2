library(Seurat)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Bseurat is the already-loaded annotated Seurat object.
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
cell_colors <- c(
  "CD8 NV"="#DFC4DE", "CD8 hobit"="#b092b2",
  "CD8 mem. 1"="#d475a8", "CD8 mem. 2"="#a8c5db",
  "CD8 mem. 3"="#7eb5d6", "CD8 eff. 1"="#569ecf",
  "CD8 eff. 2"="#3d7ab8", "CD8 IFN"="#6b429e",
  "CD4 NV"="#e38888", "CD4 mem."="#e36442",
  "Th17"="#d9382b", "Treg"="#a81c1c", "CD4 IFN"="#7a1122",
  "γδT"="#8b52a3", "MAIT"="#5c2d82",
  "MLP"="#1f8a8a", "MEP"="#229eab", "NP"="#3bc4be",
  "MDP"="#42b883", "pDCs"="#7bc993", "cDCs"="#4db86e",
  "preB"="#2e9e52", "pro/pre T"="#1b7d3a", "MP/mono"="#115c27"
)

genes <- c(
 'CD3G','TRBC2','CD8A','CD8B','CCR7','SELL','ZNF683','GZMK',
 'CD160','NKG7','PRF1','GZMB','GZMH','GNLY','FGFBP2','ADGRG1','ISG15','IFI6',
 'CD4','IL2','FAM13A','CCR10','IL21','CTSH','LMNA','FOXP3','IL2RA',
 'TRDC','TRGC1','KLRB1','SLC4A10','NCR3',
 'CD34','GATA2','CTSG','AZU1','IRF7','SCT','TGFBI','DNTT','VPREB1','ITGAM'
)
cell_order <- names(cell_colors)
x <- GetAssayData(Bseurat, assay='RNA', slot='data')
celltype <- as.character(Bseurat$Celltype)
celltype[celltype=='gdT'] <- 'γδT'
avg <- matrix(NA_real_, length(genes),length(cell_order),dimnames=list(genes,cell_order))
for (ct in cell_order) {
 ix <- which(celltype==ct)
 avg[,ct] <- log1p(Matrix::rowMeans(expm1(x[genes,ix,drop=FALSE])))
}
m <- t(scale(t(avg)))
m[!is.finite(m)] <- 0
m[m < -2] <- -2
m[m > 2] <- 2
colour_fun <- colorRamp2(c(-2,0,2),c('#326A82','#FBF9F4','#B52635'))
column_family <- factor(c(rep('CD8',8),rep('CD4',5),'γδT','MAIT',rep('HSPC',9)),levels=c('CD8','CD4','γδT','MAIT','HSPC'))
gene_family <- factor(c(rep('CD8',18),rep('CD4',9),rep('γδT',2),rep('MAIT',3),rep('HSPC',10)),levels=c('CD8','CD4','γδT','MAIT','HSPC'))
top <- HeatmapAnnotation(Celltype=cell_order,col=list(Celltype=cell_colors),show_legend=FALSE,show_annotation_name=FALSE,simple_anno_size=unit(2.5,'mm'))
hm <- Heatmap(
 m, name='Expression\n(z-score)', col=colour_fun,
 cluster_rows=FALSE, cluster_columns=FALSE,
 row_split=gene_family, cluster_row_slices=FALSE,row_title=NULL,
 column_split=column_family,cluster_column_slices=FALSE,column_title=c('CD8','CD4','','','HSPC'),
 top_annotation=top,row_gap=unit(1.5,'mm'),column_gap=unit(c(2,3,6,3),'mm'),
 column_title_gp=gpar(fontfamily='Arial',fontsize=15,fontface='bold'),
 row_names_side='left',
 row_names_gp=gpar(fontfamily='Arial',fontsize=14,fontface=ifelse(genes%in%c('FGFBP2','ADGRG1'),'bold','plain')),
 column_names_gp=gpar(fontfamily='Arial',fontsize=13),column_names_rot=90,
 rect_gp=gpar(col='#D8D8D8',lwd=.35),border=TRUE,
 width=unit(ncol(m)*5.6,'mm'),height=unit(nrow(m)*5.6,'mm'),
 heatmap_legend_param=list(title_gp=gpar(fontfamily='Arial',fontsize=15,fontface='bold'),labels_gp=gpar(fontfamily='Arial',fontsize=14),at=c(-2,0,2))
)
draw(hm, heatmap_legend_side='right')

# Native Cairo export, using a new filename:
# cairo_pdf('marker_heatmap_new.pdf',width=10.2,height=12.2,family='Arial')
# draw(hm,heatmap_legend_side='right')
# dev.off()
