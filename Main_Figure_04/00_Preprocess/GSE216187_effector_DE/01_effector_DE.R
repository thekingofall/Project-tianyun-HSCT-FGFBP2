library(Seurat)
library(Matrix)
library(data.table)
library(jsonlite)
setwd('/data5/tiany/GSE216187/GSE216187_effector_DE_heatmap_v20261029_013019')
sink('logs/01_effector_DE.log',split=TRUE)
tryCatch({
cat('START ',format(Sys.time()),'\n')
cat('ComplexHeatmap available: ',requireNamespace('ComplexHeatmap',quietly=TRUE),'\n')
cat('circlize available: ',requireNamespace('circlize',quietly=TRUE),'\n')
if (!exists('GSE216187_effector_source',inherits=FALSE))
  GSE216187_effector_source <- readRDS('/data5/tiany/GSE216187/GSE216187_seurat_scRNA_SCT_orig_Celltype_v20261027.rds')
src <- GSE216187_effector_source
cells <- colnames(src)[src$Celltype%in%c('CD8 eff. 1','CD8 eff. 2')]
md <- as.data.table(src@meta.data[cells,],keep.rownames='cell')
stopifnot(setequal(unique(md$orig),c('CR1','CR2','CR3','REL1','REL2','REL3')))
stopifnot(all(as.character(md$group)==sub('[0-9]+$','',as.character(md$orig))))
md[,orig:=as.character(orig)]
count <- GetAssayData(src,assay='RNA',slot='counts')[,cells,drop=FALSE]
linear <- count %*% Diagonal(x=10000/Matrix::colSums(count))
rownames(linear)<-rownames(count);colnames(linear)<-cells
people<-c('REL1','REL2','REL3','CR1','CR2','CR3')
cell_number<-md[,.(cells=.N),by=.(orig,group,Celltype)]
fwrite(cell_number,'data/processed/cells_by_patient_subtype.tsv',sep='\t')
fwrite(md[,.(cell,orig,patient_id,group,Celltype)],'data/processed/effector_cell_metadata.tsv',sep='\t')
print(cell_number)

compare_outcomes <- function(use_cells) {
  m<-linear[,use_cells,drop=FALSE];mm<-md[match(colnames(m),cell)]
  mean_linear<-sapply(people,function(id)Matrix::rowMeans(m[,mm$orig==id,drop=FALSE]))
  rownames(mean_linear)<-rownames(m)
  logmean<-log1p(mean_linear)
  p<-vapply(seq_len(nrow(logmean)),function(i)tryCatch(t.test(logmean[i,c('CR1','CR2','CR3')],logmean[i,c('REL1','REL2','REL3')])$p.value,error=function(e)NA_real_),numeric(1))
  de<-data.table(gene=rownames(m),mean_CR=rowMeans(mean_linear[,c('CR1','CR2','CR3')]),mean_REL=rowMeans(mean_linear[,c('REL1','REL2','REL3')]),P=p,
                pct_CR=Matrix::rowMeans(m[,mm$group=='CR',drop=FALSE]>0),pct_REL=Matrix::rowMeans(m[,mm$group=='REL',drop=FALSE]>0))
  de[,log2FC:=log2((mean_CR+1)/(mean_REL+1))]
  de[,direction:=ifelse(log2FC>0,'CR higher',ifelse(log2FC<0,'REL higher','Equal'))]
  de[,DEG:=!is.na(P)&P<.05&abs(log2FC)>=.25&pmax(pct_CR,pct_REL)>=.10]
  list(DE=de,logmean=logmean,mean_linear=mean_linear)
}
res<-compare_outcomes(cells)
fwrite(res$DE[order(P)],'data/processed/CR_vs_REL_combined_effectors_DE.tsv',sep='\t')
eff2<-compare_outcomes(md[Celltype=='CD8 eff. 2']$cell)
fwrite(eff2$DE[order(P)],'data/processed/CR_vs_REL_eff2_only_DE.tsv',sep='\t')
markers_CR<-c('KLRG1','IKZF1','GZMH','CX3CR1','ZNF683','FGFBP2','KLRD1')
markers_REL<-c('CXCR4','DUSP2','CD27')
markers<-c(markers_CR,markers_REL)
stopifnot(all(markers%in%res$DE$gene))
target<-res$DE[match(markers,gene)];target[,requested_group:=rep(c('CR','REL'),c(7,3))]
sens<-eff2$DE[match(markers,gene),.(gene,eff2_log2FC=log2FC,eff2_P=P)]
target<-merge(target,sens,by='gene',sort=FALSE)
fwrite(target,'data/processed/requested_gene_results.tsv',sep='\t');print(target)

technical <- grepl('^MT-|^RP[SL][0-9]|^TR[ABDG][VJ]|^IG[HKL][VDJ]',res$DE$gene)|res$DE$gene%in%c('XIST','TSIX','RPS4Y1','DDX3Y','EIF1AY','KDM5D','UTY','ZFY')
plot_de<-copy(res$DE[!technical & DEG]);plot_de[,magnitude:=abs(log2FC)]
setorder(plot_de,P,-magnitude,gene)
selected<-unique(c(head(plot_de[log2FC>0]$gene,100),head(plot_de[log2FC<0]$gene,100),markers))
zmat<-t(scale(t(res$logmean[selected,,drop=FALSE])))
zmat[!is.finite(zmat)]<-0
fwrite(data.table(gene=rownames(zmat),zmat),'data/plot_ready/heatmap_row_zscores.tsv',sep='\t')
fwrite(data.table(gene=rownames(res$logmean),res$logmean),'data/plot_ready/patient_log_expression.tsv',sep='\t')
gene_info<-res$DE[match(selected,gene)];gene_info[,selection:=ifelse(gene%in%markers,'Requested marker','DEG')]
fwrite(gene_info,'data/plot_ready/heatmap_gene_info.tsv',sep='\t')
sample_info<-md[,.(group=as.character(group[1]),cells=.N,eff1_cells=sum(Celltype=='CD8 eff. 1'),eff2_cells=sum(Celltype=='CD8 eff. 2')),by=orig]
fwrite(sample_info,'data/plot_ready/patient_annotation.tsv',sep='\t')
saveRDS(list(DE=res$DE,patient_expression=res$logmean,eff2_DE=eff2$DE),'data/processed/effector_analysis_cache.rds')
write_json(list(input='/data5/tiany/GSE216187/GSE216187_seurat_scRNA_SCT_orig_Celltype_v20261027.rds',cells=length(cells),individuals=6,CR_cells=sum(md$group=='CR'),REL_cells=sum(md$group=='REL'),tested_genes=sum(is.finite(res$DE$P)),CR_up_DEGs=sum(plot_de$log2FC>0),REL_up_DEGs=sum(plot_de$log2FC<0),heatmap_genes=length(selected),markers_CR=markers_CR,markers_REL=markers_REL),'reports/analysis_audit.json',pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file='logs/sessionInfo.txt')
cat('EFFECTOR_DE_READY\n')
},error=function(e){cat('EFFECTOR_DE_ERROR ',conditionMessage(e),'\n');traceback()})
sink()
