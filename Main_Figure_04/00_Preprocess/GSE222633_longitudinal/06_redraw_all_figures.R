# Run from the version directory, using only the compact exported tables.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
local({
redraw_dir <- paste0("visualization_", format(Sys.time(), "%Y%m%d_%H%M%S"))
if (dir.exists(redraw_dir)) stop("Output directory already exists")
dir.create(file.path(redraw_dir,"figures/png"), recursive=TRUE)
dir.create(file.path(redraw_dir,"figures/pdf"), recursive=TRUE)
if(!file.exists('data/plot_ready/Teff_gene_expression.csv'))stop('Set working directory to the delivered version folder')
if(Sys.info()[['sysname']]=='Linux')Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')
pal<-fread('palette_registry.tsv');cols<-setNames(pal$hex,pal$Pair)
order<-c('Pre_GCSF','Post_GCSF','D90','D180');labs<-c('Pre','Post','D90','D180')
prep<-function(d){d[,`:=`(x=match(timepoint,order),Pair=factor(Pair,levels=LETTERS[1:5]))];d}
th<-theme_classic(base_family='Arial',base_size=17)+theme(axis.title=element_text(face='bold'),axis.text=element_text(color='black'),axis.line=element_line(linewidth=.85,color='black'),axis.ticks=element_line(linewidth=.7),legend.position='right',legend.title=element_text(face='bold'),strip.background=element_blank(),strip.text=element_text(face='bold',size=17),plot.margin=margin(12,14,12,12))
axis<-scale_x_continuous(breaks=1:4,labels=labs,limits=c(.85,4.15))
# Full grids retain missing samples as NA. Segments only connect adjacent observed stages.
trajectory<-function(d,value,ylab,fac=NULL){
 bycols<-c('Pair',fac);d<-prep(copy(d));d[,value:=get(value)];grid<-CJ(Pair=factor(LETTERS[1:5],levels=LETTERS[1:5]),x=1:4)
 if(length(fac)){grid<-CJ(Pair=factor(LETTERS[1:5],levels=LETTERS[1:5]),x=1:4,facet_value=unique(d[[fac]]));setnames(grid,'facet_value',fac)}
 q<-merge(grid,d[,c(bycols,'x','value'),with=FALSE],by=c(bycols,'x'),all.x=TRUE);setorderv(q,c(bycols,'x'));q[,`:=`(xend=shift(x,type='lead'),yend=shift(value,type='lead')),by=bycols]
 seg<-q[is.finite(value)&is.finite(yend)];seg[,kind:=ifelse(x==2,'Donor to recipient','Within person')]
 p<-ggplot(q,aes(x,value,color=Pair))+geom_segment(data=seg,aes(xend=xend,yend=yend,linetype=kind),linewidth=.85,show.legend=c(linetype=FALSE))+geom_segment(data=seg[Pair=='D'],aes(xend=xend,yend=yend,linetype=kind),linewidth=1.35,show.legend=FALSE)+geom_point(aes(shape=Pair),size=3.3,stroke=.6,na.rm=TRUE)+scale_color_manual(values=cols,labels=paste('Pair',LETTERS[1:5]))+scale_shape_manual(values=c(16,17,15,18,8),labels=paste('Pair',LETTERS[1:5]))+scale_linetype_manual(values=c('Donor to recipient'='dashed','Within person'='solid'))+axis+labs(x=NULL,y=ylab,color=NULL,shape=NULL)+th
 if(length(fac))p<-p+facet_wrap(as.formula(paste('~',fac)),nrow=1,scales='free_y')
 p
}
manifest<-list()
savefig<-function(p,name,w=8,h=5){ggsave(file.path(redraw_dir,paste0('figures/png/',name,'.png')),p,width=w,height=h,dpi=300,bg='white',type='cairo');ggsave(file.path(redraw_dir,paste0('figures/pdf/',name,'.pdf')),p,width=w,height=h,device=cairo_pdf,family='Arial');manifest[[length(manifest)+1]]<<-data.table(figure=name,width=w,height=h,png=file.path(redraw_dir,paste0('figures/png/',name,'.png')),pdf=file.path(redraw_dir,paste0('figures/pdf/',name,'.pdf')))}
ex<-fread('data/plot_ready/Teff_gene_expression.csv')
for(g in c('FGFBP2','ADGRG1')){p<-trajectory(ex[gene==g],'mean_log_expression',paste0(g,' expression\n(mean log-normalized RNA)'));savefig(p,paste0('Fig',ifelse(g=='FGFBP2',1,2),'_',g,'_Teff_trajectory'));p2<-trajectory(ex[gene==g],'positive_percent',paste0(g,'-positive Teff (%)'))+scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100));savefig(p2,paste0('FigS',ifelse(g=='FGFBP2',1,2),'_',g,'_positive_trajectory'))}
po<-fread('data/plot_ready/Teff_positive_percent.csv');savefig(trajectory(po,'Double_positive','Double-positive Teff (%)')+scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100)),'FigS3_double_positive_trajectory')
co<-fread('data/plot_ready/sample_composition.csv');savefig(trajectory(co,'Teff_percent','Teff among CD8 T cells (%)')+scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100)),'FigS4_Teff_composition')
sc<-fread('data/plot_ready/Teff_functional_scores.csv');sl<-melt(sc,id.vars=c('sample','Pair','participant_id','role','timepoint'),variable.name='Program',value.name='score');sl[,Program:=factor(Program,levels=c('Cytotoxicity','Exhaustion','Naive_memory'),labels=c('Cytotoxicity','Exhaustion','Naive / memory'))];savefig(trajectory(sl,'score','Mean signature expression','Program'),'FigS5_independent_functional_programs',12.8,4.8)
ps<-fread('data/plot_ready/donor_Teff_recipient_detection.csv');ps[,timepoint:=factor(timepoint,levels=c('Pre_GCSF','Post_GCSF'),labels=c('Pre','Post'))];ps[,Expression_group:=factor(Expression_group,levels=c('Other Teff','Double positive'))];p<-ggplot(ps,aes(Expression_group,recipient_detected_percent,group=Pair,color=factor(Pair,levels=LETTERS[1:5])))+geom_line(linewidth=.8)+geom_point(aes(shape=factor(Pair,levels=LETTERS[1:5])),size=3)+facet_wrap(~timepoint,nrow=1)+scale_color_manual(values=cols,labels=paste('Pair',LETTERS[1:5]))+scale_shape_manual(values=c(16,17,15,18,8),labels=paste('Pair',LETTERS[1:5]))+scale_y_continuous(limits=c(0,100),breaks=c(0,25,50,75,100))+labs(x=NULL,y='Donor cells with clonotypes\ndetected in recipients (%)',color=NULL,shape=NULL)+th+theme(axis.text.x=element_text(size=14));savefig(p,'FigS6_TCR_recipient_detection',9.5,5)
fwrite(rbindlist(manifest),file.path(redraw_dir,'figure_manifest.tsv'),sep='\t');cat('DONE_PLOT',nrow(rbindlist(manifest)),'figures\n')
})

