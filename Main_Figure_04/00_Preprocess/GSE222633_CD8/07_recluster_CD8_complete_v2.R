suppressPackageStartupMessages({library(Seurat);library(harmony);library(data.table);library(Matrix);library(UCell);library(jsonlite)})
setwd('/data5/tiany/GSE222633/GSE222633_CD8_annotation_v20261029_132636')
stopifnot(!file.exists('logs/07_recluster_CD8_complete_v2.log'))
logfile<-file('logs/07_recluster_CD8_complete_v2.log','wt');sink(logfile,split=TRUE);sink(logfile,type='message')
tryCatch({
set.seed(222633)
cat('RECLUSTER_START ',format(Sys.time()),'\n');flush(logfile)
parent_CD8<-get('GSE222633',envir=.GlobalEnv)
CD8_source_cells<-colnames(parent_CD8)[as.character(parent_CD8$celltype_broad)=='CD8 T']
stopifnot(length(CD8_source_cells)==53538L)
FCD8_work<-subset(parent_CD8,cells=CD8_source_cells)
FCD8_work$source_seurat_clusters<-as.character(FCD8_work$seurat_clusters)
FCD8_work@reductions<-list();FCD8_work@graphs<-list();FCD8_work@neighbors<-list()
DefaultAssay(FCD8_work)<-'RNA'
FCD8_work<-NormalizeData(FCD8_work,normalization.method='LogNormalize',scale.factor=10000,verbose=FALSE)
FCD8_work<-FindVariableFeatures(FCD8_work,selection.method='vst',nfeatures=3000,verbose=FALSE)
features<-VariableFeatures(FCD8_work)
features<-features[!grepl('^MT-|^RP[SL][0-9]|^TR[ABDG][VJ]|^IG[HKL][VDJ]',features)]
features<-setdiff(features,c('XIST','TSIX','RPS4Y1','DDX3Y','EIF1AY','KDM5D','UTY','ZFY'))
VariableFeatures(FCD8_work)<-features
fwrite(data.table(gene=features),'data/processed/CD8_variable_features.tsv',sep='\t')
FCD8_work<-ScaleData(FCD8_work,features=features,verbose=FALSE)
FCD8_work<-RunPCA(FCD8_work,features=features,npcs=50,seed.use=222633,verbose=FALSE)
cat('PCA_READY ',format(Sys.time()),' features=',length(features),'\n');flush(logfile)
harmony_method<-getS3method('RunHarmony','Seurat',envir=asNamespace('harmony'))
FCD8_work<-harmony_method(object=FCD8_work,group.by.vars='pair_id',reduction='pca',dims.use=1:30,
 max.iter.harmony=20,project.dim=FALSE,assay.use='RNA',verbose=TRUE)
FCD8_work<-FindNeighbors(FCD8_work,reduction='harmony',dims=1:30,graph.name=c('CD8_nn','CD8_snn'),verbose=FALSE)
FCD8_work<-FindClusters(FCD8_work,graph.name='CD8_snn',resolution=.5,algorithm=1,random.seed=222633,verbose=FALSE)
FCD8_work$CD8_cluster<-as.character(Idents(FCD8_work))
FCD8_work<-RunUMAP(FCD8_work,reduction='harmony',dims=1:30,n.neighbors=30,min.dist=.3,seed.use=222633,verbose=FALSE)
cat('CLUSTER_UMAP_READY ',format(Sys.time()),'\n');print(table(FCD8_work$CD8_cluster));flush(logfile)
md<-as.data.table(FCD8_work@meta.data,keep.rownames='cell')
emb<-as.data.table(Embeddings(FCD8_work,'umap'),keep.rownames='cell');setnames(emb,c('cell','UMAP_1','UMAP_2'))
fwrite(merge(md[,.(cell,CD8_cluster,sample,pair_id,participant_id,role,timepoint,TCR_paired_TRA_TRB,TCR_unambiguous_TRA_TRB)],emb,by='cell'),
 'data/plot_ready/CD8_UMAP_preannotation.tsv',sep='\t')
fwrite(md[,.(cells=.N,strict_TCR_percent=100*mean(TCR_unambiguous_TRA_TRB%in%TRUE),paired_TCR_percent=100*mean(TCR_paired_TRA_TRB%in%TRUE)),by=CD8_cluster],
 'data/processed/CD8_cluster_sizes_TCR.tsv',sep='\t')
fwrite(md[,.(cells=.N),by=.(CD8_cluster,sample,pair_id,role,timepoint)],'data/processed/CD8_cluster_sample_counts.tsv',sep='\t')

sets<-fread('resources/annotation_gene_sets.tsv');programs<-lapply(split(sets$gene,sets$program),intersect,rownames(FCD8_work))
x<-GetAssayData(FCD8_work,assay='RNA',slot='data')
scores<-as.matrix(ScoreSignatures_UCell(x,features=programs,ncores=2,chunk.size=500,maxRank=1500))
stopifnot(identical(rownames(scores),colnames(FCD8_work)))
colnames(scores)<-paste0('score_',sub('_UCell$','',colnames(scores)))
FCD8_work<-AddMetaData(FCD8_work,metadata=as.data.frame(scores))
fwrite(as.data.table(FCD8_work@meta.data[,c('CD8_cluster',colnames(scores))])[,lapply(.SD,mean),by=CD8_cluster,.SDcols=colnames(scores)],
 'data/processed/CD8_cluster_program_scores.tsv',sep='\t')
groups<-sort(unique(FCD8_work$CD8_cluster));cl<-match(FCD8_work$CD8_cluster,groups);ns<-tabulate(cl,nbins=length(groups))
W<-sparseMatrix(i=seq_along(cl),j=cl,x=1/ns[cl],dims=c(length(cl),length(groups)))
panel<-unique(c(sets$gene,'CD8A','CD8B','CD4','CD74','HLA-DRA','HLA-DRB1','CD38','TNFRSF9','NR4A1','NR4A2','FOS','JUN','BCL2','GZMA','GZMM','KLRD1','GZMK','S100A4','S100A6','MS4A1','CD79A','LYZ','LST1','PPBP','HBB'))
panel<-intersect(panel,rownames(x))
lin<-x[panel,,drop=FALSE];lin@x<-expm1(lin@x)
avg<-log1p(as.matrix(lin%*%W));colnames(avg)<-groups
pct<-100*as.matrix((lin>0)%*%W);colnames(pct)<-groups
long_avg<-as.data.table(as.table(avg));setnames(long_avg,c('gene','CD8_cluster','expression'))
long_pct<-as.data.table(as.table(pct));setnames(long_pct,c('gene','CD8_cluster','percent_detected'))
fwrite(merge(long_avg,long_pct,by=c('gene','CD8_cluster')),'data/plot_ready/CD8_cluster_marker_panel.tsv',sep='\t')
adt<-GetAssayData(FCD8_work,assay='ADT',slot='data');am<-as.matrix(adt%*%W);colnames(am)<-groups
alt<-as.data.table(as.table(am));setnames(alt,c('protein','CD8_cluster','ADT_CLR'))
fwrite(alt,'data/processed/CD8_cluster_ADT_profiles.tsv',sep='\t')
cat('ANNOTATION_PROFILES_READY ',format(Sys.time()),'\n');flush(logfile)

Idents(FCD8_work)<-FCD8_work$CD8_cluster
CD8_markers_work<-FindAllMarkers(FCD8_work,assay='RNA',slot='data',only.pos=TRUE,min.pct=.1,
 logfc.threshold=.25,test.use='wilcox',max.cells.per.ident=1500,random.seed=222633,return.thresh=.05,verbose=TRUE)
stopifnot(nrow(CD8_markers_work)>0)
fwrite(as.data.table(CD8_markers_work),'data/processed/CD8_FindAllMarkers_unannotated.csv')
fwrite(as.data.table(CD8_markers_work)[order(cluster,-avg_log2FC),head(.SD,25),by=cluster],
 'data/processed/CD8_top25_markers_by_cluster.tsv',sep='\t')
capture.output(sessionInfo(),file='logs/sessionInfo.txt')
write_json(list(status='markers_ready',cells=ncol(FCD8_work),RNA_genes=nrow(FCD8_work[['RNA']]),clusters=length(groups),
 HVGs=length(features),dims=30,seed=222633,cluster_counts=as.list(table(FCD8_work$CD8_cluster)),marker_rows=nrow(CD8_markers_work)),
 'reports/recluster_audit.json',pretty=TRUE,auto_unbox=TRUE)
cat('CD8_RECLUSTER_MARKERS_READY ',format(Sys.time()),'\n');flush(logfile)
},error=function(e){cat('RECLUSTER_ERROR ',conditionMessage(e),'\n');traceback();flush(logfile)})
sink(type='message');sink();close(logfile)
