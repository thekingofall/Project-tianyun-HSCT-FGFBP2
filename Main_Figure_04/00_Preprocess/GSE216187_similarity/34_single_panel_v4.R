library(data.table)
library(ggplot2)
path <- '/data5/tiany/GSE216187/S2_similarity_v20261028_215325'
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
d <- fread(file.path(path,'data/plot_ready/similarity_rankings.tsv'))
d <- d[method=='Centered HVG, markers excluded'][order(rank)]
stopifnot(nrow(d)==8L)
d[,celltype:=factor(celltype,levels=rev(celltype))]
palette <- fread(file.path(path,'palette_registry.tsv'))
cell_colors <- setNames(palette$hex,palette$label)

p <- ggplot(d,aes(rho,celltype,colour=celltype))+
  geom_segment(aes(x=0,xend=rho,yend=celltype),linewidth=1,alpha=.65)+
  geom_point(size=4.5)+
  geom_text(aes(label=sprintf('%.3f',rho)),nudge_y=.25,
            colour='black',family='Arial',size=4.9)+
  scale_colour_manual(values=cell_colors,guide='none')+
  scale_x_continuous(expand=expansion(mult=c(.10,.12)))+
  labs(x='Spearman correlation with CD8_S2_Teff',y=NULL)+
  theme_classic(base_family='Arial',base_size=19)+
  theme(axis.line=element_line(linewidth=.8,colour='black'),
        axis.ticks=element_line(linewidth=.7,colour='black'),
        axis.text=element_text(size=16,colour='black'),
        axis.title.x=element_text(size=19),
        plot.margin=margin(14,16,12,12))

png <- file.path(path,'figures/png/S2_01_centered_single_panel_v4.png')
pdf <- file.path(path,'figures/pdf/S2_01_centered_single_panel_v4.pdf')
table <- file.path(path,'data/plot_ready/S2_centered_single_panel_v4.tsv')
stopifnot(!file.exists(png),!file.exists(pdf),!file.exists(table))
fwrite(d,table,sep='\t')
ggsave(png,p,width=9,height=7.2,dpi=300,device=grDevices::png,type='cairo',bg='white')
ggsave(pdf,p,width=9,height=7.2,device=grDevices::cairo_pdf,family='Arial',bg='white')
cat('SINGLE_PANEL_V4_READY\n')
