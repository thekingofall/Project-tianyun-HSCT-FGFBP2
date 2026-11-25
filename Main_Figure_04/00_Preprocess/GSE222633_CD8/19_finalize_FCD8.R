suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite)})
local({
root<-'/data5/tiany/GSE222633/GSE222633_CD8_annotation_v20261029_132636';setwd(root)
stopifnot(!file.exists('logs/19_finalize_FCD8.log'))
logfile<-file('logs/19_finalize_FCD8.log','wt');sink(logfile,split=TRUE);sink(logfile,type='message')
on.exit({sink(type='message');sink();close(logfile)},add=TRUE)
cat('FINALIZE_START ',format(Sys.time()),'\n');flush(logfile)
FCD8<-get('FCD8_work',envir=.GlobalEnv)
parent<-get('GSE222633',envir=.GlobalEnv)
ann<-fread('resources/CD8_annotation_map.tsv');ann[,cluster:=as.character(cluster)]
stopifnot(setequal(unique(FCD8$CD8_cluster),ann$cluster),ncol(FCD8)==53538L)
ii<-match(as.character(FCD8$CD8_cluster),ann$cluster);stopifnot(!anyNA(ii))
FCD8$Celltype<-ann$Celltype[ii]
FCD8$Celltype_fine<-ann$Celltype_fine[ii]
FCD8$Celltype2<-FCD8$Celltype_fine
FCD8$annotation_confidence<-ann$confidence[ii]
FCD8$annotation_review<-as.logical(ann$review[ii])
FCD8$annotation_lineage<-ann$lineage[ii]
FCD8$CD8_annotation_evidence<-ann$evidence[ii]
FCD8$is_conventional_CD8<-FCD8$annotation_lineage=='Conventional CD8 T'
FCD8$CD8_high_confidence<-FCD8$is_conventional_CD8 & !FCD8$annotation_review
FCD8$seurat_clusters<-factor(FCD8$CD8_cluster)
type_order<-c('CD8 Naive','CD8 Tem','CD8 Teff','CD8 Cycling','CD8 Activated','MAIT-like T','gdT-like T','CD8 T (low RNA)','CD4-like T (review)','NK/T-like (review)')
FCD8$Celltype<-factor(FCD8$Celltype,levels=type_order)
FCD8$timepoint<-factor(FCD8$timepoint,levels=c('Pre_GCSF','Post_GCSF','D90','D180'))
clinical<-fread('resources/pair_clinical_metadata.tsv');jj<-match(as.character(FCD8$pair_id),as.character(clinical$pair_id));stopifnot(!anyNA(jj))
FCD8$Pair<-factor(clinical$Pair[jj],levels=LETTERS[1:5])
FCD8$recipient_diagnosis<-clinical$recipient_diagnosis[jj]
FCD8$participant_age<-ifelse(FCD8$role=='Donor',clinical$donor_age[jj],clinical$recipient_age[jj])
FCD8$participant_sex<-ifelse(FCD8$role=='Donor',clinical$donor_sex[jj],clinical$recipient_sex[jj])
FCD8$pair_recipient_aGVHD_grade<-clinical$aGVHD_overall_grade[jj]
FCD8$pair_recipient_cGVHD_report<-clinical$cGVHD[jj]
FCD8$pair_recipient_relapse_report<-clinical$relapse[jj]
FCD8@misc$CD8_reanalysis<-list(source='/data5/tiany/GSE222633/GSE222633_seurat_RNA_TCR_harmony_v20261026.rds',
  source_selection='celltype_broad == CD8 T',source_cells=53538L,seed=222633,HVGs=length(VariableFeatures(FCD8)),
  PCA_dimensions=50L,Harmony_group='pair_id',Harmony_dimensions=30L,main_resolution=.5,subcluster_resolution=.4,
  annotation_map=as.data.frame(ann),annotation_genes=as.data.frame(fread('resources/annotation_gene_sets.tsv')),
  classical_reference='Manuscript Naive/Tem/Teff programs; Mathioudaki et al Blood2024 Fig1; Obermayer et al FrontImmunol2023',
  annotation_note='Tcm, Tex, IFN and tissue residence were not forced when no independent marker-supported population was resolved. Innate-like and review categories retained.',
  TCR_note='All original TCR and clone metadata retained; inherited frequency columns keep their original source definitions.',
  marker_method='FindAllMarkers Wilcoxon/limma; only.pos TRUE; min.pct0.1; logfc.threshold0.25; max.cells.per.ident1500; random.seed222633')
DefaultAssay(FCD8)<-'RNA';Idents(FCD8)<-FCD8$Celltype
stopifnot(nrow(FCD8[['RNA']])==33538L,all(FCD8$celltype_broad=='CD8 T'),!anyNA(FCD8$Celltype))
tcrcols<-grep('^TCR|^CT|^clone|^clonal|^Frequency$',colnames(parent@meta.data),value=TRUE)
source_meta<-parent@meta.data[colnames(FCD8),tcrcols,drop=FALSE]
target_meta<-FCD8@meta.data[,tcrcols,drop=FALSE]
tcr_ok<-all(vapply(tcrcols,function(k)identical(as.character(source_meta[[k]]),as.character(target_meta[[k]])),logical(1)))
stopifnot(tcr_ok,identical(colnames(FCD8),CD8_source_cells))
raw_src<-GetAssayData(parent,assay='RNA',slot='counts')[,colnames(FCD8),drop=FALSE]
raw_dst<-GetAssayData(FCD8,assay='RNA',slot='counts')
counts_ok<-identical(raw_src@p,raw_dst@p)&&identical(raw_src@i,raw_dst@i)&&identical(raw_src@x,raw_dst@x)
stopifnot(counts_ok);validObject(FCD8)
md<-as.data.table(FCD8@meta.data,keep.rownames='cell')
emb<-as.data.table(Embeddings(FCD8,'umap'),keep.rownames='cell');setnames(emb,c('cell','UMAP_1','UMAP_2'))
umap<-merge(md[,.(cell,Celltype=as.character(Celltype),Celltype_fine,Celltype2,CD8_cluster,CD8_cluster_initial,annotation_confidence,annotation_review,sample,Pair,role,timepoint)],emb,by='cell')
fwrite(umap,'data/plot_ready/FCD8_UMAP_annotated.tsv',sep='\t')
fwrite(md[,.(cells=.N),by=.(Celltype,Celltype_fine,CD8_cluster,annotation_confidence,annotation_review)],'data/processed/FCD8_annotation_counts.tsv',sep='\t')
fwrite(md[,.(cells=.N),by=.(sample,Pair,participant_id,role,timepoint,Celltype)],'data/processed/FCD8_annotation_by_sample.tsv',sep='\t')
fwrite(md[,.(cells=.N),by=Celltype],'data/plot_ready/FCD8_celltype_counts.tsv',sep='\t')
print(table(FCD8$Celltype));flush(logfile)
genes<-c('CD3D','CD8A','CD4','CCR7','TCF7','LEF1','IL7R','GZMK','CCL5','FGFBP2','GZMB','GZMH','ZNF683','MKI67','TOP2A','PDCD1','LAG3','TOX','KLRB1','SLC4A10','TRDC','TRGC2','FCER1G')
genes<-intersect(genes,rownames(FCD8));x<-GetAssayData(FCD8,assay='RNA',slot='data')[genes,,drop=FALSE];x@x<-expm1(x@x)
g<-as.character(FCD8$Celltype);groups<-type_order[type_order%in%g];ci<-match(g,groups);ns<-tabulate(ci,nbins=length(groups))
W<-sparseMatrix(i=seq_along(ci),j=ci,x=1/ns[ci],dims=c(length(ci),length(groups)))
av<-log1p(as.matrix(x%*%W));colnames(av)<-groups
pc<-100*as.matrix((x>0)%*%W);colnames(pc)<-groups
sc<-t(scale(t(av)));sc[!is.finite(sc)]<-0
a<-as.data.table(as.table(av));setnames(a,c('gene','Celltype','mean_log_expression'))
b<-as.data.table(as.table(pc));setnames(b,c('gene','Celltype','percent_detected'))
c<-as.data.table(as.table(sc));setnames(c,c('gene','Celltype','scaled_expression'))
fwrite(Reduce(function(x,y)merge(x,y,by=c('gene','Celltype')),list(a,b,c)),'data/plot_ready/FCD8_annotation_marker_dotplot.tsv',sep='\t')
markers<-as.data.table(get('CD8_markers_final_work',envir=.GlobalEnv));markers[,cluster:=as.character(cluster)]
markers<-merge(markers,ann[,.(cluster,Celltype,Celltype_fine,confidence,review)],by='cluster',sort=FALSE)
stopifnot(nrow(markers)>0)
out_rds<-'/data5/tiany/GSE222633/FCD8.rds';out_markers<-'/data5/tiany/GSE222633/CD8diff.csv'
if(file.exists(out_rds)||file.exists(out_markers)){
 out_rds<-file.path(root,'FCD8.rds');out_markers<-file.path(root,'CD8diff.csv')
 cat('Existing parent targets retained; using version-directory targets.\n')
}
stopifnot(!file.exists(out_rds),!file.exists(out_markers))
fwrite(markers,out_markers)
assign('FCD8',FCD8,envir=.GlobalEnv)
cat('SAVING_RDS ',out_rds,' ',format(Sys.time()),'\n');flush(logfile)
saveRDS(FCD8,out_rds,compress='gzip')
cat('RDS_SAVED ',file.info(out_rds)$size,' bytes ',format(Sys.time()),'\n');flush(logfile)
check<-readRDS(out_rds)
stopifnot(ncol(check)==53538L,nrow(check[['RNA']])==33538L,!anyNA(check$Celltype),
 all(c('RNA','ADT','HTO')%in%Assays(check)),all(c('pca','harmony','umap')%in%Reductions(check)),
 ncol(Embeddings(check,'harmony'))==30L,ncol(Embeddings(check,'umap'))==2L,
 identical(as.character(check$Celltype),as.character(FCD8$Celltype)))
validObject(check)
audit<-list(status='PASS',RDS=out_rds,markers_csv=out_markers,cells=ncol(check),RNA_genes=nrow(check[['RNA']]),
 assays=Assays(check),reductions=Reductions(check),initial_clusters=length(unique(check$CD8_cluster_initial)),
 refined_clusters=length(unique(check$CD8_cluster)),celltype_counts=as.list(table(check$Celltype)),
 TCR_columns_preserved=length(tcrcols),TCR_values_unchanged=tcr_ok,RNA_counts_unchanged=counts_ok,
 strict_TRA_TRB_cells=sum(check$TCR_unambiguous_TRA_TRB%in%TRUE),review_cells=sum(check$annotation_review),
 marker_rows=nrow(markers),RDS_bytes=file.info(out_rds)$size,fresh_readRDS_valid=TRUE)
write_json(audit,'reports/FCD8_delivery_readback.json',pretty=TRUE,auto_unbox=TRUE)
cat('FCD8_DELIVERY_VALIDATED ',format(Sys.time()),'\n');flush(logfile)
})
