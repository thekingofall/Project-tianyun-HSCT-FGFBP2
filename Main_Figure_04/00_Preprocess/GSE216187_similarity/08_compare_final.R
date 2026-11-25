library(Seurat)
library(Matrix)
library(data.table)
library(UCell)
library(jsonlite)
setwd('/data5/tiany/GSE216187/S2_similarity_v20261028_215325')
sink('logs/08_compare_final.log',split=TRUE)
tryCatch({
set.seed(216187)
cat('COMPARE_START ',format(Sys.time()),'\n')
ext_cells <- colnames(External)[grepl('^CD8 ',External$Celltype)]
ref_cells <- colnames(HSCT)[grepl('^S[123]_CD8_',HSCT$Celltype2)]
E <- CreateSeuratObject(GetAssayData(External,assay='RNA',slot='counts')[,ext_cells],meta.data=External@meta.data[ext_cells,])
R <- CreateSeuratObject(GetAssayData(HSCT,assay='RNA',slot='counts')[,ref_cells],meta.data=HSCT@meta.data[ref_cells,])
E <- FindVariableFeatures(NormalizeData(E,verbose=FALSE),nfeatures=2000,verbose=FALSE)
R <- FindVariableFeatures(NormalizeData(R,verbose=FALSE),nfeatures=2000,verbose=FALSE)
common <- intersect(rownames(E),rownames(R))
nuisance <- unique(c(grep('^MT-|^RP[SL][0-9]|^TR[ABDG][VJ]|^IG[HKL][VDJ]',common,value=TRUE),
                    c('XIST','TSIX','RPS4Y1','DDX3Y','EIF1AY','KDM5D','UTY','ZFY')))
genes_all <- setdiff(intersect(VariableFeatures(E),VariableFeatures(R)),nuisance)
genes <- setdiff(genes_all,c('FGFBP2','ADGRG1'))
genes_union <- setdiff(intersect(union(VariableFeatures(E),VariableFeatures(R)),common),c(nuisance,'FGFBP2','ADGRG1'))
stopifnot(length(genes)>200)
fwrite(data.table(gene=common),'data/processed/shared_genes.tsv',sep='\t')
fwrite(data.table(gene=genes),'data/processed/primary_features.tsv',sep='\t')
fwrite(data.table(gene=genes_union),'data/processed/union_features.tsv',sep='\t')

em <- data.table(cell=colnames(E),celltype=E$Celltype,sample=E$sample_id,subject=E$patient_id,alias=E$orig,group=E$group)
rm <- data.table(cell=colnames(R),celltype=as.character(R$Celltype2),sample=as.character(R$orig.ident),subject=sub('[ACD]$','',as.character(R$orig.ident)),group=as.character(R$group))
fwrite(em,'data/processed/external_CD8_metadata.tsv',sep='\t');fwrite(rm,'data/processed/HSCT_CD8_metadata.tsv',sep='\t')
e_order <- c('CD8 NV','CD8 hobit','CD8 mem. 1','CD8 mem. 2','CD8 mem. 3','CD8 eff. 1','CD8 eff. 2','CD8 IFN')
r_order <- c('S1_CD8_Naive','S1_CD8_Tcm','S1_CD8_Tem','S1_CD8_Teff','S1_CD8_Tex','S2_CD8_Naive','S2_CD8_Tcm','S2_CD8_Tem','S2_CD8_Teff','S2_CD8_Tex','S3_CD8_Naive','S3_CD8_Tcm','S3_CD8_Tem','S3_CD8_Teff','S3_CD8_Tex')
ref_types <- r_order[r_order%in%rm[,.N,by=celltype][N>=100]$celltype]

# Small sample/subtype means keep cell counts from becoming sample weights.
make_profiles <- function(ob,md,types) {
  md <- copy(md);md[,key:=paste(subject,sample,celltype,sep='|')]
  units <- md[,.(n=.N,subject=subject[1],sample=sample[1],celltype=celltype[1]),by=key][n>=10 & celltype%in%types]
  x <- GetAssayData(ob,assay='RNA',slot='data')[common,,drop=FALSE];x@x<-expm1(x@x)
  ind <- match(md$key,units$key);ok<-which(!is.na(ind))
  W <- sparseMatrix(i=ok,j=ind[ok],x=1/units$n[ind[ok]],dims=c(ncol(x),nrow(units)))
  sample_mean <- log1p(as.matrix(x%*%W));rownames(sample_mean)<-common;colnames(sample_mean)<-units$key
  units[,subject_key:=paste(subject,celltype,sep='|')]
  su <- units[,.(subject=subject[1],celltype=celltype[1],samples=.N,cells=sum(n)),by=subject_key]
  W2 <- sparseMatrix(i=seq_len(nrow(units)),j=match(units$subject_key,su$subject_key),x=1/su$samples[match(units$subject_key,su$subject_key)],dims=c(nrow(units),nrow(su)))
  subject_mean <- as.matrix(sample_mean%*%W2);rownames(subject_mean)<-common;colnames(subject_mean)<-su$subject_key
  balance <- sapply(types,function(ct)rowMeans(subject_mean[,su$celltype==ct,drop=FALSE]));rownames(balance)<-common
  pooled <- sapply(types,function(ct)log1p(Matrix::rowMeans(x[,md$celltype==ct,drop=FALSE])));rownames(pooled)<-common
  list(balanced=balance,pooled=pooled,subject=subject_mean,subject_meta=su,sample=sample_mean,sample_meta=units)
}
EP <- make_profiles(E,em,e_order);RP <- make_profiles(R,rm,ref_types)
saveRDS(list(external=EP,reference=RP,genes=genes,genes_all=genes_all,genes_union=genes_union), 'data/processed/profile_cache.rds')
fwrite(data.table(gene=common,EP$balanced),'data/plot_ready/external_balanced_profiles.tsv',sep='\t')
fwrite(data.table(gene=common,RP$balanced),'data/plot_ready/reference_balanced_profiles.tsv',sep='\t')
fwrite(EP$subject_meta,'data/processed/external_subject_units.tsv',sep='\t');fwrite(RP$subject_meta,'data/processed/HSCT_subject_units.tsv',sep='\t')
center <- function(m)sweep(m,1,rowMeans(m),'-')
RC <- center(RP$balanced);EC <- center(EP$balanced)
cross_cor <- cor(RC[genes,],EC[genes,],method='spearman')
corr <- as.data.table(as.table(cross_cor));setnames(corr,c('reference','celltype','rho'))
fwrite(corr,'data/plot_ready/centered_similarity_heatmap.tsv',sep='\t')
target <- 'S2_CD8_Teff'
metrics <- list(
  'Centered HVG, markers excluded'=cor(RC[genes,target],EC[genes,],method='spearman')[1,],
  'Centered HVG, markers included'=cor(RC[genes_all,target],EC[genes_all,],method='spearman')[1,],
  'Centered union HVG'=cor(RC[genes_union,target],EC[genes_union,],method='spearman')[1,],
  'Uncentered HVG'=cor(RP$balanced[genes,target],EP$balanced[genes,],method='spearman')[1,],
  'Centered pooled cells'=cor(center(RP$pooled)[genes,target],center(EP$pooled)[genes,],method='spearman')[1,])
ranking <- rbindlist(lapply(names(metrics),function(nm)data.table(method=nm,celltype=names(metrics[[nm]]),rho=as.numeric(metrics[[nm]]))))
ranking[,rank:=frank(-rho,ties.method='min'),by=method]
fwrite(ranking,'data/plot_ready/similarity_rankings.tsv',sep='\t')
print(ranking[method=='Centered HVG, markers excluded'][order(rank)])

# Perturb participants, keeping the gene set fixed; these are stability checks.
stability<-list();idx<-0L
for (rd in unique(RP$subject_meta$subject)) for (ed in unique(EP$subject_meta$subject)) {
  rb <- sapply(ref_types,function(ct){ii<-RP$subject_meta$celltype==ct & RP$subject_meta$subject!=rd;if(!any(ii))return(rep(NA_real_,length(common)));rowMeans(RP$subject[,ii,drop=FALSE])})
  eb <- sapply(e_order,function(ct){ii<-EP$subject_meta$celltype==ct & EP$subject_meta$subject!=ed;if(!any(ii))return(rep(NA_real_,length(common)));rowMeans(EP$subject[,ii,drop=FALSE])})
  rownames(rb)<-rownames(eb)<-common
  if(anyNA(rb)||anyNA(eb))next
  rho<-as.numeric(cor(center(rb)[genes,target],center(eb)[genes,],method='spearman'))
  idx<-idx+1L;stability[[idx]]<-data.table(omit_HSCT=rd,omit_external=ed,celltype=e_order,rho=rho,rank=frank(-rho,ties.method='min'))
}
stability<-rbindlist(stability);fwrite(stability,'data/processed/leave_subject_out_runs.tsv',sep='\t')
stab <- stability[,.(runs=.N,top1=sum(rank==1),top1_percent=100*mean(rank==1),median_rho=median(rho),min_rho=min(rho),max_rho=max(rho)),by=celltype]
fwrite(stab,'data/plot_ready/ranking_stability.tsv',sep='\t')

# Relearn a complete S2 state signature from the manuscript object itself.
delta<-RP$balanced[,target]-rowMeans(RP$balanced[,setdiff(ref_types,target),drop=FALSE])
detect <- Matrix::rowMeans(GetAssayData(R,slot='counts')[common,R$Celltype2==target,drop=FALSE]>0)
eligible<-setdiff(common,c(nuisance,'FGFBP2','ADGRG1'))
up<-head(names(sort(delta[eligible][detect[eligible]>=.10 & delta[eligible]>0],decreasing=TRUE)),100)
down<-head(names(sort(delta[eligible][delta[eligible]<0],decreasing=FALSE)),100)
stopifnot(length(up)>=20,length(down)>=20)
signature<-data.table(gene=c(up,down),direction=rep(c('S2 up','S2 down'),c(length(up),length(down))),contrast=delta[c(up,down)])
fwrite(signature,'data/processed/S2_signature_genes.tsv',sep='\t')
cat('SIGNATURE_READY ',length(up),length(down),'\n')
us <- ScoreSignatures_UCell(GetAssayData(E,slot='data')[common,],features=list(S2_up=up,S2_down=down),maxRank=1500,ncores=2,chunk.size=500)
em[,S2_signature:=us[cell,'S2_up_UCell']-us[cell,'S2_down_UCell']]
fwrite(em,'data/processed/external_CD8_signature_cells.tsv',sep='\t')
sig_subject<-em[,.(score=mean(S2_signature),cells=.N,alias=alias[1],group=group[1]),by=.(subject,celltype)][cells>=10]
fwrite(sig_subject,'data/plot_ready/S2_signature_subject_scores.tsv',sep='\t')
sig_summary<-sig_subject[,.(signature_mean=mean(score),signature_sd=sd(score),subjects=.N),by=celltype]
fwrite(sig_summary,'data/plot_ready/S2_signature_summary.tsv',sep='\t')
heatgenes<-unique(c('FGFBP2','ADGRG1',head(up,20),head(down,10)))
hm <- cbind(RP$balanced[heatgenes,c('S2_CD8_Teff','S3_CD8_Teff','S3_CD8_Tex')],EP$balanced[heatgenes,])
z<-t(scale(t(hm)));z[!is.finite(z)]<-0
ht<-as.data.table(as.table(z));setnames(ht,c('gene','state','z'))
fwrite(ht,'data/plot_ready/signature_expression_heatmap.tsv',sep='\t')
fwrite(data.table(order=seq_along(heatgenes),gene=heatgenes),'data/plot_ready/heatmap_gene_order.tsv',sep='\t')
audit<-list(external_cells=ncol(E),HSCT_CD8_cells=ncol(R),target_cells=sum(R$Celltype2==target),shared_genes=length(common),primary_genes=length(genes),union_genes=length(genes_union),signature_up=length(up),signature_down=length(down),reference_states=ref_types,external_states=e_order,stability_runs=idx,external_input='/data5/tiany/GSE216187/GSE216187_seurat_scRNA_SCT_orig_Celltype_v20261027.rds',HSCT_input='/data5/tiany/HSCT/Result/Tcell.rds')
write_json(audit,'reports/analysis_audit.json',pretty=TRUE,auto_unbox=TRUE)
capture.output(sessionInfo(),file='logs/sessionInfo.txt')
cat('COMPARE_READY ',format(Sys.time()),'\n')
},error=function(e){cat('COMPARE_ERROR ',conditionMessage(e),'\n');traceback()})
sink()
