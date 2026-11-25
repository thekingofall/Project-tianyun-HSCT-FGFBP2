suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite);library(digest)})

finalize_gse222633<-function(){
 output<-'/data5/tiany/GSE222633/GSE222633_seurat_RNA_TCR_harmony_v20261026.rds'
 stopifnot(exists('GSE222633',envir=.GlobalEnv),!file.exists(output),!file.exists('reports/Seurat_delivery_verified.json'))
 obj<-get('GSE222633',envir=.GlobalEnv)
 download<-read_json('data/raw/download_manifest.json',simplifyVector=FALSE)
 stopifnot(identical(download$status,'complete'),length(download$files)==35L)
 if(file.exists('data/processed/cluster_annotation.tsv')){
  annotations<-fread('data/processed/cluster_annotation.tsv')
  stopifnot(!anyDuplicated(annotations$cluster),all(as.character(obj$seurat_clusters)%in%annotations$cluster))
  ix<-match(as.character(obj$seurat_clusters),annotations$cluster)
  obj$celltype_broad<-annotations$celltype_broad[ix];obj$celltype_annotation_evidence<-annotations$evidence[ix]
  obj@misc$celltype_annotation<-as.data.frame(annotations)
 }
 obj@misc$provenance<-list(dataset='GSE222633',paper='10.3389/fimmu.2023.1114368',source='https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE222633',author_code='https://github.com/bihealth/obermayer_et_al_tcell_persistence',source_files=download$files,user_selected_integration='Harmony by donor-recipient pair; preserve timepoints',processing='RNA/ADT/HTO quality control, DoubletFinder, LogNormalize, variable features excluding VDJ, PCA, Harmony, neighbor graph, multiresolution Louvain clustering, UMAP, native scRepertoire receptor integration',new_Seurat_file=output,software=read_json('logs/R_environment.json',simplifyVector=TRUE),time=format(Sys.time(),'%FT%T%z'))
 qc<-as.data.table(obj@misc$RNA_processing$QC_summary)
 stopifnot(ncol(obj)==sum(qc$retained_cells),identical(colnames(obj),rownames(obj@meta.data)),!anyDuplicated(colnames(obj)),all(c('pca','harmony','umap')%in%Reductions(obj)),nlevels(obj$pair_id)==5)
 for(r in Reductions(obj))stopifnot(identical(rownames(Embeddings(obj,r)),colnames(obj)),all(is.finite(Embeddings(obj,r))))
 before<-object_signature(obj)
 summary<-list(status='saving',output=output,cells=ncol(obj),genes=nrow(obj),samples=length(unique(obj$sample)),pairs=nlevels(obj$pair_id),assays=Assays(obj),reductions=Reductions(obj),graphs=names(obj@graphs),clusters=nlevels(Idents(obj)),timepoints=as.list(table(obj$timepoint)),TCR_any=sum(obj$TCR_any_valid),TCR_paired=sum(obj$TCR=='TCR+'),TCR_unambiguous_paired=sum(obj$TCR_unambiguous_TRA_TRB),TCR_multichain=sum(obj$TCR_multichain),persisting_clonotypes=sum(obj@misc$TCR$paired_clone_summary$persisting),persisting_cells=sum(obj$TCR_persisting),account=Sys.info()[['user']],pid=Sys.getpid())
 stage('save_final_RNA_TCR_Seurat')
 con<-gzfile(output,'wb',compression=1);tryCatch(saveRDS(obj,con,version=3),finally=close(con))
 write_json(summary,'logs/Seurat_saved.json',auto_unbox=TRUE,pretty=TRUE)
 stage('independently_reload_and_verify_saved_object')
 check<-readRDS(output)
 stopifnot(identical(before,object_signature(check)))
 methods::validObject(check)
 sum_coverage<-sum(check@misc$TCR$coverage$paired_TRA_TRB)
 stopifnot(sum_coverage==sum(check$TCR_paired_TRA_TRB))
 summary$status<-'complete';summary$bytes<-unname(file.info(output)$size);summary$sha256<-digest(output,file=TRUE,algo='sha256');summary$fresh_readback<-'PASS';summary$all_assays_counts_normalized_scaled_values<-'PASS';summary$metadata_and_cell_order<-'PASS';summary$reductions_graphs_and_identities<-'PASS';summary$source_file_count<-length(download$files);summary$time<-format(Sys.time(),'%FT%T%z')
 write_json(summary,'reports/Seurat_delivery_verified.json',auto_unbox=TRUE,pretty=TRUE)
 writeLines(paste(summary$sha256,output),'reports/Seurat_delivery.sha256')
 fwrite(as.data.table(check@meta.data,keep.rownames='cell_id'),'data/processed/final_cell_metadata.tsv.gz',sep='\t')
 writeLines(capture.output(sessionInfo()),'logs/final_sessionInfo.txt')
 assign('GSE222633',check,envir=.GlobalEnv)
 if(exists('GSE222633_QC_samples',envir=.GlobalEnv))rm(list='GSE222633_QC_samples',envir=.GlobalEnv)
 gc()
 stage('FINAL_SEURAT_DELIVERED_AND_VERIFIED')
 print(summary);invisible(summary)
}

object_signature<-function(o){
 assays<-setNames(lapply(Assays(o),function(a){x<-o[[a]];list(counts=digest(x@counts,algo='sha256'),data=digest(x@data,algo='sha256'),scale_data=digest(x@scale.data,algo='sha256'),features=digest(rownames(x),algo='sha256'))}),Assays(o))
 list(dim=dim(o),cells=digest(colnames(o),algo='sha256'),genes=digest(rownames(o),algo='sha256'),assays=assays,metadata=lapply(o@meta.data,digest,algo='sha256'),metadata_row_names=digest(rownames(o@meta.data),algo='sha256'),reductions=lapply(Reductions(o),function(r)digest(o[[r]],algo='sha256')),graphs=lapply(o@graphs,digest,algo='sha256'),identities=digest(Idents(o),algo='sha256'),misc=digest(o@misc,algo='sha256'),variable_features=digest(VariableFeatures(o),algo='sha256'))
}
