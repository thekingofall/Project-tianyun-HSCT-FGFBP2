suppressPackageStartupMessages({library(Seurat);library(jsonlite);library(data.table);library(digest)})
verification_source<-readLines("/data5/tiany/GSE222633/GSE222633_v20261026_222500/code/15_memory_bounded_verification.R",warn=FALSE);eval(parse(text=verification_source[-1L]))

save_gse216187<-function(){
 output<-'/data5/tiany/GSE216187/GSE216187_seurat_scRNA_SCT_v20261027.rds'
 stopifnot(exists('GSE216187',envir=.GlobalEnv),!file.exists(output),!file.exists('reports/Seurat_delivery_verified.json'))
 obj<-get('GSE216187',envir=.GlobalEnv)
 obj$tissue<-"Bone marrow";obj$timepoint<-"D100";obj$diagnosis<-"AML"
 if(file.exists('data/processed/cluster_annotation.tsv')){
  a<-fread('data/processed/cluster_annotation.tsv');ix<-match(as.character(obj$seurat_clusters),as.character(a$cluster));stopifnot(!anyNA(ix))
  obj$celltype_broad<-a$celltype_broad[ix];obj$celltype_annotation_evidence<-a$evidence[ix];obj@misc$cluster_annotation<-as.data.frame(a)
 }
 obj$tissue<-'Bone marrow';obj$timepoint<-'D100'
 source_files<-read_json('data/raw/download_manifest.json',simplifyVector=FALSE)
 stopifnot(source_files$status=='complete',length(source_files$files)==8L,all(obj$souporcell_status=='singlet'),all(obj$nFeature_RNA>=200),all(obj$percent.mt<=15),!anyNA(obj$patient_id))
 stopifnot(all(obj$group[obj$GSM=='GSM6660161']=='REL'),all(obj$group[obj$GSM=='GSM6660162']=='CR'))
 obj@misc$provenance<-list(dataset='GSE216187',paper='10.1182/blood.2023021815',source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE216187',author_code='https://github.com/ZauggGroup/scBM_alloSCT',GEO_metadata_correction='2026-03-08: GSM6660161=REL; GSM6660162=CR; data contents unchanged',source_files=source_files$files,modality='scRNA-seq; SCT and RNA_uncorrected are derived representations of the same RNA data',software=setNames(lapply(c('Seurat','Matrix','sctransform','SoupX','org.Hs.eg.db'),function(p)as.character(packageVersion(p))),c('Seurat','Matrix','sctransform','SoupX','org.Hs.eg.db')),time=format(Sys.time(),'%FT%T%z'))
 baseline<-list(core=object_signature(obj),SCT_models=digest(obj[['SCT']]@SCTModel.list,algo='sha256'))
 qc<-as.data.table(obj@misc$processing$QC_summary)
 stopifnot(ncol(obj)==sum(qc$retained),identical(colnames(obj),rownames(obj@meta.data)),!anyDuplicated(colnames(obj)))
 stage('save_final_scRNA_Seurat')
 con<-gzfile(output,'wb',compression=1);tryCatch(saveRDS(obj,con,version=3),finally=close(con))
 stage('reload_and_verify_saved_scRNA_object')
 check<-readRDS(output)
 stopifnot(identical(baseline,list(core=object_signature(check),SCT_models=digest(check[['SCT']]@SCTModel.list,algo='sha256'))))
 methods::validObject(check)
 result<-list(status='complete',output=output,cells=ncol(check),RNA_genes=nrow(check[['RNA']]),assays=Assays(check),reductions=Reductions(check),graphs=names(check@graphs),clusters=nlevels(Idents(check)),patients=nlevels(check$patient_id),group_counts=as.list(table(check$group)),source_files=8,bytes=unname(file.info(output)$size),sha256=digest(output,file=TRUE,algo='sha256'),fresh_readback='PASS',all_count_normalized_scaled_matrices='PASS',SCT_models='PASS',metadata_reductions_graphs_identities='PASS',GEO_corrected_group_labels='PASS',account=Sys.info()[['user']],pid=Sys.getpid(),time=format(Sys.time(),'%FT%T%z'))
 write_json(result,'reports/Seurat_delivery_verified.json',auto_unbox=TRUE,pretty=TRUE);writeLines(paste(result$sha256,output),'reports/Seurat_delivery.sha256')
 fwrite(as.data.table(check@meta.data,keep.rownames='cell_id'),'data/processed/final_cell_metadata.tsv.gz',sep='\t')
 writeLines(capture.output(sessionInfo()),'logs/final_sessionInfo.txt')
 assign('GSE216187',check,envir=.GlobalEnv)
 if(exists('GSE216187_QC_samples',envir=.GlobalEnv))rm(list='GSE216187_QC_samples',envir=.GlobalEnv)
 gc();stage('FINAL_SC_RNA_OBJECT_VERIFIED_COMPLETE');print(result);invisible(result)
}
