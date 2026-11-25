suppressPackageStartupMessages({library(data.table);library(ggplot2);library(ggrepel)})
setwd('/data5/tiany/GSE222633/GSE222633_CD8_annotation_v20261029_132636')
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
dir.create('figures/png',recursive=TRUE,showWarnings=FALSE);dir.create('figures/pdf',recursive=TRUE,showWarnings=FALSE)
palette<-fread('palette_registry.tsv');cols<-setNames(palette$hex,palette$Celltype);type_order<-palette$Celltype
labels_clean<-function(x)gsub(' \\(review\\)','',x)
theme_cd8<-theme_classic(base_family='Arial',base_size=19)+theme(axis.line=element_line(linewidth=.8,colour='black'),
 axis.ticks=element_line(linewidth=.7,colour='black'),axis.text=element_text(size=15,colour='black'),
 legend.title=element_blank(),legend.text=element_text(size=14),plot.margin=margin(10,12,10,10))
save_pair<-function(p,id,w,h,data){
 a<-paste0('figures/png/',id,'.png');b<-paste0('figures/pdf/',id,'.pdf');stopifnot(!file.exists(a),!file.exists(b))
 ggsave(a,p,width=w,height=h,dpi=300,device=grDevices::png,type='cairo',bg='white')
 ggsave(b,p,width=w,height=h,device=grDevices::cairo_pdf,family='Arial',bg='white')
 data.table(figure_id=id,PNG=a,PDF=b,data=data,width=w,height=h)
}
d<-fread('data/plot_ready/FCD8_UMAP_annotated.tsv');d[,Celltype:=factor(Celltype,levels=type_order)]
set.seed(222633);d<-d[sample(.N)]
lab<-d[,{
 mx<-median(UMAP_1);my<-median(UMAP_2);i<-which.min((UMAP_1-mx)^2+(UMAP_2-my)^2)
 .(UMAP_1=UMAP_1[i],UMAP_2=UMAP_2[i],cells=.N)
},by=Celltype]
lab[,label:=labels_clean(as.character(Celltype))]
p1<-ggplot(d,aes(UMAP_1,UMAP_2,colour=Celltype))+geom_point(size=.25,alpha=.8,shape=16)+
 geom_label_repel(data=lab[cells>=500],aes(label=label),family='Arial',size=4.7,colour='black',fill='white',
  label.size=.25,box.padding=.45,point.padding=.1,min.segment.length=0,seed=222633,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=cols,labels=labels_clean)+labs(x='UMAP 1',y='UMAP 2')+theme_cd8+
 guides(colour=guide_legend(override.aes=list(size=3,alpha=1)))
m1<-save_pair(p1,'FCD8_UMAP_Celltype',11,8,'data/plot_ready/FCD8_UMAP_annotated.tsv')
clab<-d[,{
 mx<-median(UMAP_1);my<-median(UMAP_2);i<-which.min((UMAP_1-mx)^2+(UMAP_2-my)^2)
 .(UMAP_1=UMAP_1[i],UMAP_2=UMAP_2[i])
},by=CD8_cluster]
p2<-ggplot(d,aes(UMAP_1,UMAP_2,colour=Celltype))+geom_point(size=.25,alpha=.8,shape=16)+
 geom_label_repel(data=clab,aes(label=CD8_cluster),family='Arial',size=4.1,colour='black',fill='white',
  label.size=.2,seed=222633,max.overlaps=Inf,show.legend=FALSE)+
 scale_colour_manual(values=cols,labels=labels_clean)+labs(x='UMAP 1',y='UMAP 2')+theme_cd8+
 guides(colour=guide_legend(override.aes=list(size=3,alpha=1)))
m2<-save_pair(p2,'FCD8_UMAP_clusters',11,8,'data/plot_ready/FCD8_UMAP_annotated.tsv')
dot<-fread('data/plot_ready/FCD8_annotation_marker_dotplot.tsv')
gene_order<-c('CD3D','CD8A','CD4','CCR7','TCF7','LEF1','IL7R','GZMK','CCL5','FGFBP2','GZMB','GZMH','ZNF683','MKI67','TOP2A','PDCD1','LAG3','TOX','KLRB1','SLC4A10','TRDC','TRGC2','FCER1G')
dot[,`:=`(gene=factor(gene,levels=gene_order),Celltype=factor(Celltype,levels=rev(type_order)))]
p3<-ggplot(dot,aes(gene,Celltype,size=percent_detected,colour=scaled_expression))+geom_point()+
 scale_size_area(max_size=5,limits=c(0,100),breaks=c(25,50,75,100),name='Detected (%)')+
 scale_colour_gradient2(low='#2166AC',mid='#F7F7F7',high='#B2182B',midpoint=0,limits=c(-2,2),oob=scales::squish,name='Scaled expression')+
 scale_y_discrete(labels=labels_clean)+labs(x=NULL,y=NULL)+theme_cd8+
 theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=13),axis.text.y=element_text(size=14),legend.position='bottom')+
 guides(colour=guide_colourbar(barwidth=grid::unit(45,'mm'),barheight=grid::unit(4,'mm'),title.position='top'))
m3<-save_pair(p3,'FCD8_annotation_marker_DotPlot',13.2,7.6,'data/plot_ready/FCD8_annotation_marker_dotplot.tsv')
fwrite(rbindlist(list(m1,m2,m3)),'figure_manifest.tsv',sep='\t')
cat('FCD8_FIGURES_READY\n')
