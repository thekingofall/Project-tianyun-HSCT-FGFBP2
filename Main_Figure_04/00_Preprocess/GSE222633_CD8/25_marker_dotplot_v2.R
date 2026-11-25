suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setwd('/data5/tiany/GSE222633/GSE222633_CD8_annotation_v20261029_132636')
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
dot<-fread('data/plot_ready/FCD8_annotation_marker_dotplot.tsv')
types<-fread('palette_registry.tsv')$Celltype
genes<-c('CD3D','CD8A','CD4','CCR7','TCF7','LEF1','IL7R','GZMK','CCL5','FGFBP2','GZMB','GZMH','ZNF683','MKI67','TOP2A','PDCD1','LAG3','TOX','KLRB1','SLC4A10','TRDC','TRGC2','FCER1G')
dot[,`:=`(gene=factor(gene,levels=genes),Celltype=factor(Celltype,levels=rev(types)))]
p<-ggplot(dot,aes(gene,Celltype,size=percent_detected,colour=scaled_expression))+geom_point()+
 scale_size_area(max_size=5,limits=c(0,100),breaks=c(25,50,75,100),name='Detected (%)')+
 scale_colour_gradient2(low='#2166AC',mid='#F7F7F7',high='#B2182B',midpoint=0,limits=c(-2,2),oob=scales::squish,name='Scaled expression')+
 scale_y_discrete(labels=function(x)gsub(' \\(review\\)','',x))+labs(x=NULL,y=NULL)+
 theme_classic(base_family='Arial',base_size=19)+
 theme(axis.line=element_line(linewidth=.8,colour='black'),axis.ticks=element_line(linewidth=.7,colour='black'),
  axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=13,colour='black'),axis.text.y=element_text(size=14,colour='black'),
  legend.position='bottom',legend.title=element_text(size=14,family='Arial'),legend.text=element_text(size=14),plot.margin=margin(10,12,10,10))+
 guides(size=guide_legend(title.position='top',title.hjust=.5),
 colour=guide_colourbar(barwidth=grid::unit(45,'mm'),barheight=grid::unit(4,'mm'),title.position='top',title.hjust=.5))
a<-'figures/png/FCD8_annotation_marker_DotPlot_v2.png';b<-'figures/pdf/FCD8_annotation_marker_DotPlot_v2.pdf'
stopifnot(!file.exists(a),!file.exists(b))
ggsave(a,p,width=13.2,height=7.6,dpi=300,device=grDevices::png,type='cairo',bg='white')
ggsave(b,p,width=13.2,height=7.6,device=grDevices::cairo_pdf,family='Arial',bg='white')
cat('MARKER_DOTPLOT_V2_READY\n')
