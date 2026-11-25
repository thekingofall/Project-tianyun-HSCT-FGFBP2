suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ComplexHeatmap);library(circlize);library(grid);library(jsonlite)})
cfg<-read_json('resources/figure_runtime.json',simplifyVector=TRUE)
Sys.setenv(FONTCONFIG_FILE=cfg$fontconfig);Sys.unsetenv('MAO_ARIAL_FONT_DIR')
for(path in c('figures/png','figures/pdf'))dir.create(path,recursive=TRUE,showWarnings=FALSE)
cairo_pdf(tempfile(fileext='.pdf'),family='Arial',width=10,height=12)
group_colours<-c(CR='#3B6C80',REL='#B36350')
expression_colours<-circlize::colorRamp2(c(-2,0,2),c('#326A82','#FBF9F4','#B52635'))
patient_order<-c('CR1','CR2','CR3','REL1','REL2','REL3')
cell_palette<-fread('resources/celltype_palette.tsv');ct_order<-cell_palette$Celltype
colours<-setNames(cell_palette$hex,cell_palette$Celltype)
family<-setNames(c(rep('CD8',8),rep('CD4',5),'γδT','MAIT',rep('HSPC',9)),ct_order)
marker<-fread('resources/marker_gene_order.tsv')
comparisons<-fread('data/processed/between_group_statistics.tsv')
paired_tests<-fread('data/processed/paired_state_statistics.tsv')
manifest<-list()
p_label<-function(p)ifelse(is.na(p),'P = NA',ifelse(p<.001,'P < 0.001',paste0('P = ',formatC(p,format='g',digits=2))))

save_heatmap<-function(hm,id,w,h,data_path){
 png_path<-paste0('figures/png/',id,'.png');pdf_path<-paste0('figures/pdf/',id,'.pdf')
 stopifnot(!file.exists(png_path),!file.exists(pdf_path))
 grDevices::png(png_path,width=w,height=h,units='in',res=300,type='cairo',bg='white')
 ComplexHeatmap::draw(hm,heatmap_legend_side='right',annotation_legend_side='right',merge_legends=TRUE,padding=unit(c(5,5,5,5),'mm'));dev.off()
 grDevices::cairo_pdf(pdf_path,width=w,height=h,family='Arial',bg='white')
 ComplexHeatmap::draw(hm,heatmap_legend_side='right',annotation_legend_side='right',merge_legends=TRUE,padding=unit(c(5,5,5,5),'mm'));dev.off()
 manifest[[id]]<<-data.table(figure_id=id,png=png_path,pdf=pdf_path,width=w,height=h,dpi=300,data=data_path,renderer='code/11_render_marker_heatmap_final.R',status='Rendered')
}
save_plot<-function(p,id,w,h,data_path){
 png_path<-paste0('figures/png/',id,'.png');pdf_path<-paste0('figures/pdf/',id,'.pdf')
 stopifnot(!file.exists(png_path),!file.exists(pdf_path))
 ggsave(png_path,p,width=w,height=h,units='in',dpi=300,device=grDevices::png,type='cairo',bg='white')
 ggsave(pdf_path,p,width=w,height=h,units='in',device=grDevices::cairo_pdf,family='Arial',bg='white')
 manifest[[id]]<<-data.table(figure_id=id,png=png_path,pdf=pdf_path,width=w,height=h,dpi=300,data=data_path,renderer='code/11_render_marker_heatmap_final.R',status='Rendered')
}
theme_validation<-function(){theme_classic(base_family='Arial',base_size=19)+theme(axis.line=element_line(colour='black',linewidth=.8),axis.ticks=element_line(colour='black',linewidth=.7),axis.text=element_text(colour='black',size=16),axis.title=element_text(size=18),strip.background=element_rect(fill='white',colour='black',linewidth=.8),strip.text=element_text(size=18,face='bold'),legend.position='bottom',legend.title=element_blank(),legend.text=element_text(size=15),plot.margin=margin(10,12,8,10),panel.spacing=unit(10,'mm'))}

# 1. Marker-by-celltype expression heatmap with the requested ordering.
d<-fread('data/plot_ready/marker_celltype_heatmap.tsv')
wide<-dcast(d,gene~Celltype,value.var='plot_zscore');m<-as.matrix(wide[match(marker$gene,wide$gene),..ct_order]);rownames(m)<-marker$gene
stopifnot(!anyNA(m),identical(colnames(m),ct_order))
top<-HeatmapAnnotation(Celltype=ct_order,col=list(Celltype=colours),show_legend=FALSE,show_annotation_name=FALSE,simple_anno_size=unit(2.5,'mm'))
hm<-Heatmap(m,name='Expression\n(z-score)',col=expression_colours,cluster_rows=FALSE,cluster_columns=FALSE,
 row_split=factor(marker$gene_family,levels=unique(marker$gene_family)),cluster_row_slices=FALSE,row_title=NULL,row_gap=unit(1.5,'mm'),
 column_split=factor(family[ct_order],levels=c('CD8','CD4','γδT','MAIT','HSPC')),cluster_column_slices=FALSE,column_gap=unit(c(2,3,6,3),'mm'),
 column_title=c('CD8','CD4','','','HSPC'),column_title_gp=gpar(fontfamily='Arial',fontsize=15,fontface='bold'),top_annotation=top,
 row_names_side='left',row_names_gp=gpar(fontfamily='Arial',fontsize=14,fontface=ifelse(rownames(m)%in%c('FGFBP2','ADGRG1'),'bold','plain')),
 column_names_gp=gpar(fontfamily='Arial',fontsize=13),column_names_rot=90,
 rect_gp=gpar(col='#D8D8D8',lwd=.35),border=TRUE,width=unit(ncol(m)*5.6,'mm'),height=unit(nrow(m)*5.6,'mm'),
 heatmap_legend_param=list(title_gp=gpar(fontfamily='Arial',fontsize=15,fontface='bold'),labels_gp=gpar(fontfamily='Arial',fontsize=14),at=c(-2,0,2)))
save_heatmap(hm,'Validation_01_marker_celltype_heatmap_v3',10.2,12.2,'data/plot_ready/marker_celltype_heatmap.tsv')

fwrite(rbindlist(manifest),'figure_manifest_marker_v3.tsv',sep='\t')
graphics.off();cat('MARKER_HEATMAP_V3_READY\n')
