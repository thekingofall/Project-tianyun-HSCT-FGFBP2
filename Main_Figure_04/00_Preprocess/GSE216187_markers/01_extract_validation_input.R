suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite)})
local({
 root<-'/data5/tiany/GSE216187/GSE216187_external_validation_v20261027_222329'
 oldwd<-getwd();on.exit(setwd(oldwd));setwd(root)
 logfile<-file('logs/01_extract_input.log','wt');sink(logfile,split=TRUE);on.exit({sink();close(logfile)},add=TRUE)
 stopifnot(!file.exists('data/processed/validation_input.rds'))
 obj<-get('GSE216187_annotated',envir=.GlobalEnv)
 stopifnot(ncol(obj)==8501L)
 x<-GetAssayData(obj,assay='RNA',slot='data');counts<-GetAssayData(obj,assay='RNA',slot='counts')
 md<-as.data.table(obj@meta.data,keep.rownames='cell_id')
 md[,`:=`(Celltype=as.character(Celltype),orig=as.character(orig),group=as.character(group))]
 md[Celltype=='gdT',Celltype:='γδT']
 signatures<-split(fread('resources/state_gene_sets.tsv')$gene,fread('resources/state_gene_sets.tsv')$program)
 stopifnot(!any(c('FGFBP2','ADGRG1')%in%unlist(signatures)))
 signatures<-lapply(signatures,intersect,rownames(x));stopifnot(all(lengths(signatures)>=5L))
 cd8<-grepl('^CD8',md$Celltype)
 score<-as.matrix(UCell::ScoreSignatures_UCell(x[,cd8,drop=FALSE],features=signatures,ncores=2))
 colnames(score)<-sub('_UCell$','',colnames(score))
 stopifnot(identical(rownames(score),md$cell_id[cd8]))
 for(program in names(signatures)){md[[program]]<-NA_real_;set(md,which(cd8),program,score[,program])}
 clinical_genes<-c('FGFBP2','ADGRG1','GZMB','PRF1','GZMH','GNLY','NKG7','KLRD1','CX3CR1','CD27','GZMK','CCR7','TCF7','PDCD1','LAG3','HAVCR2','CTLA4','TIGIT','TOX')
 genes<-unique(c(fread('resources/marker_gene_order.tsv')$gene,clinical_genes,unlist(signatures)))
 missing<-setdiff(genes,rownames(x));genes<-intersect(genes,rownames(x));stopifnot(!length(setdiff(c('FGFBP2','ADGRG1'),genes)))
 unc<-GetAssayData(obj,assay='RNA_uncorrected',slot='counts')[c('FGFBP2','ADGRG1'),,drop=FALSE]
 input<-list(expression=x[genes,,drop=FALSE],counts=counts[genes,,drop=FALSE],uncorrected_targets=unc,metadata=as.data.frame(md),clinical_genes=clinical_genes,signatures=signatures,missing_genes=missing,source='/data5/tiany/GSE216187/GSE216187_seurat_scRNA_SCT_orig_Celltype_v20261027.rds',score_method=paste0('UCell ',packageVersion('UCell'),'; maxRank=1500; scores calculated once on all CD8 cells; targets excluded'))
 saveRDS(input,'data/processed/validation_input.rds')
 packages<-c('Seurat','Matrix','UCell','ComplexHeatmap','circlize','ggplot2','patchwork')
 versions<-setNames(lapply(packages,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else NA_character_),packages)
 write_json(list(status='complete',cells=ncol(x),CD8_cells=sum(cd8),genes=genes,missing_genes=missing,signatures=signatures,packages=versions,groups=as.list(table(md$group)),patient_counts=as.list(table(md$orig)),score_method=input$score_method,time=format(Sys.time(),'%FT%T%z')),'reports/input_audit.json',auto_unbox=TRUE,pretty=TRUE)
 fwrite(md,'data/processed/validation_metadata.tsv.gz',sep='\t')
 cat('VALIDATION_INPUT_READY\n')
})
