library(data.table)
library(ggplot2)
path <- '/data5/tiany/GSE216187/S2_similarity_v20261028_215325'
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
d <- fread(file.path(path,'data/plot_ready/full15_similarity_heatmap_v3.tsv'))
row_order <- c('S1_Naive','S1_Tcm','S1_Tem','S1_Teff','S1_Tex',
               'S2_Naive','S2_Tcm','S2_Tem','S2_Teff','S2_Tex',
               'S3_Naive','S3_Tcm','S3_Tem','S3_Teff','S3_Tex')
col_order <- c('CD8 eff. 2','CD8 eff. 1','CD8 mem. 3','CD8 mem. 2',
               'CD8 mem. 1','CD8 IFN','CD8 hobit','CD8 NV')
d[,Subset:=factor(gsub('_CD8_','_',reference),levels=rev(row_order))]
d[,celltype:=factor(celltype,levels=col_order)]
highlight <- d[Subset=='S2_Teff' & celltype%in%c('CD8 eff. 1','CD8 eff. 2')]
stopifnot(nrow(d)==120L,nrow(highlight)==2L,!anyNA(d$rho))
p <- ggplot(d,aes(celltype,Subset,fill=rho))+
  geom_tile(colour='white',linewidth=.4)+
  geom_tile(data=highlight,fill=NA,colour='black',linewidth=1.4)+
  scale_fill_gradient2(low='#2166AC',mid='#F7F7F7',high='#B2182B',
                       midpoint=0,limits=c(-1,1),name='Spearman rho')+
  labs(x='GSE216187 CD8 subtype',y='CD8 T cell subsets')+
  theme_classic(base_family='Arial',base_size=19)+
  theme(axis.text=element_text(size=16,colour='black'),
        axis.text.x=element_text(angle=40,hjust=1,vjust=1),
        axis.ticks=element_blank(),axis.line=element_blank(),
        legend.position='bottom',legend.title=element_text(size=16),
        legend.text=element_text(size=15),plot.margin=margin(12,14,12,12))+
  guides(fill=guide_colourbar(barwidth=grid::unit(65,'mm'),
                             barheight=grid::unit(4,'mm'),title.position='top'))

p
