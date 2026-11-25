library(Seurat)
library(Matrix)
library(data.table)
setwd('/data5/tiany/GSE216187/S2_similarity_v20261028_215325')
sink('logs/28_prepare_full15_v3.log',split=TRUE)
tryCatch({
cache <- readRDS('data/processed/profile_cache.rds')
old <- fread('data/plot_ready/centered_similarity_heatmap.tsv')
full_order <- as.vector(t(outer(paste0('S',1:3,'_CD8_'),c('Naive','Tcm','Tem','Teff','Tex'),paste0)))
missing_types <- setdiff(full_order,unique(old$reference))
if (!exists('HSCT_full',inherits=FALSE)) HSCT_full <- readRDS('/data5/tiany/HSCT/Result/Tcell.rds')
md <- as.data.table(HSCT_full@meta.data,keep.rownames='cell')
md[,celltype:=as.character(Celltype2)]
md[,subject:=sub('[ACD]$','',as.character(orig.ident))]
md[,sample:=as.character(orig.ident)]
counts <- GetAssayData(HSCT_full,assay='RNA',slot='counts')
common <- rownames(cache$reference$balanced)
features <- cache$genes
ref_center <- rowMeans(cache$reference$balanced)
ext_centered <- sweep(cache$external$balanced,1,rowMeans(cache$external$balanced),'-')
new_rows<-list();unit_rows<-list()
for (ct in missing_types) {
  cm <- copy(md[celltype==ct]);stopifnot(nrow(cm)>0)
  cc <- counts[,cm$cell,drop=FALSE]
  linear <- cc %*% Diagonal(x=10000/Matrix::colSums(cc))
  rownames(linear)<-rownames(cc);colnames(linear)<-colnames(cc)
  # Include every observed cell for the four requested rare states.
  cm[,key:=paste(subject,sample,sep='|')]
  units<-cm[,.(subject=subject[1],sample=sample[1],cells=.N),by=key]
  sample_profiles<-sapply(units$key,function(k)log1p(Matrix::rowMeans(linear[common,cm$key==k,drop=FALSE])))
  rownames(sample_profiles)<-common
  subject_profiles<-sapply(unique(units$subject),function(id)rowMeans(sample_profiles[,units$subject==id,drop=FALSE]))
  rownames(subject_profiles)<-common
  balanced<-rowMeans(subject_profiles)
  # Project new rows against the frozen original 11-state reference center.
  relative<-balanced-ref_center
  rho<-as.numeric(cor(relative[features],ext_centered[features,],method='spearman'))
  new_rows[[ct]]<-data.table(reference=ct,celltype=colnames(ext_centered),rho=rho)
  units[,reference:=ct];unit_rows[[ct]]<-units
}
full<-rbindlist(c(list(old),new_rows))
stopifnot(nrow(full)==15L*8L,setequal(unique(full$reference),full_order))
fwrite(full,'data/plot_ready/full15_similarity_heatmap_v3.tsv',sep='\t')
fwrite(rbindlist(unit_rows),'data/processed/full15_added_sample_units_v3.tsv',sep='\t')
numbers<-md[celltype%in%full_order,.(cells=.N,subjects=uniqueN(subject)),by=celltype]
numbers[,order:=match(celltype,full_order)];setorder(numbers,order)
fwrite(numbers,'data/processed/full15_cell_counts_v3.tsv',sep='\t')
print(numbers)
cat('FULL15_DATA_READY\n')
},error=function(e)cat('FULL15_DATA_ERROR ',conditionMessage(e),'\n'))
sink()
