library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
Sys.setenv(FONTCONFIG_FILE='/data5/tiany/GSE295359/Fig2ABC_v20261025_135600/resources/fontconfig.xml')

em_types <- c('CD8 mem. 1','CD8 mem. 2','CD8 mem. 3','CD8 eff. 1','CD8 eff. 2')
patients <- c('CR1','CR2','CR3','REL1','REL2','REL3')
targets <- c('FGFBP2','ADGRG1')
x <- GetAssayData(Bseurat,assay='RNA',slot='data')
meta <- Bseurat@meta.data
parts <- list()
for(patient in patients) {
 ix <- which(meta$orig==patient & meta$Celltype%in%em_types)
 means <- log1p(Matrix::rowMeans(expm1(x[targets,ix,drop=FALSE])))
 parts[[patient]] <- data.table(orig=patient,group=unique(as.character(meta$group[ix])),
                               gene=targets,log1p_mean=as.numeric(means))
}
patient_expression <- rbindlist(parts)
patient_expression[,`:=`(group=factor(group,levels=c('CR','REL')),gene=factor(gene,levels=targets))]
tests <- patient_expression[,.(P=t.test(log1p_mean[group=='CR'],log1p_mean[group=='REL'])$p.value),by=gene]
tests[,label:=ifelse(P<.001,'P < 0.001',paste0('P = ',formatC(P,format='g',digits=2)))]
annotation <- patient_expression[,.(y=max(log1p_mean)*1.15+.03),by=gene]
annotation <- merge(annotation,tests,by='gene')

p <- ggplot(patient_expression,aes(group,log1p_mean,colour=group))+
 geom_point(size=4.2,position=position_jitter(width=.06,height=0,seed=216187))+
 stat_summary(fun=mean,fun.min=mean,fun.max=mean,geom='crossbar',width=.42,linewidth=.7,colour='black')+
 geom_text(data=annotation,aes(x=1.5,y=y,label=label),inherit.aes=FALSE,family='Arial',size=5.7,colour='black')+
 facet_wrap(~gene,nrow=1,scales='free_y')+
 scale_colour_manual(values=c(CR='#3B6C80',REL='#B36350'))+
 scale_y_continuous(expand=expansion(mult=c(.08,.14)))+
 labs(x=NULL,y='Mean expression\n(log1p normalized counts)')+
 theme_classic(base_family='Arial',base_size=19)+
 theme(legend.position='none',axis.text=element_text(colour='black',size=16),
       axis.line=element_line(linewidth=.8),
       strip.background=element_rect(fill='white',colour='black',linewidth=.8),
       strip.text=element_text(face='bold',size=18))
p
tests
