suppressPackageStartupMessages({library(Seurat);library(Matrix);library(SoupX);library(data.table);library(jsonlite);library(AnnotationDbi);library(org.Hs.eg.db)})
setDTthreads(2);set.seed(216187);options(warn=1,future.globals.maxSize=8*1024^3)
future::plan('sequential')

process_gse216187<-function(){
 plan<-read_json('resources/download_plan.json',simplifyVector=TRUE)
 groups<-data.table(GSM=c('GSM6660161','GSM6660162'),group=c('REL','CR'),legacy=c('CR','REL'))
 objects<-list();raw_qc<-list();summaries<-list();soup_profiles<-list();donor_tables<-list();feature_map<-NULL
 for(i in seq_len(nrow(groups))){
  gsm<-groups$GSM[i];group<-groups$group[i];stage(paste0('read_raw_matrix ',gsm,' ',group))
  paths<-plan[plan$GSM==gsm,]
  get_path<-function(pattern)paths$path[grepl(pattern,paths$filename)]
  features<-fread(get_path('features'),header=FALSE);barcodes<-scan(gzfile(get_path('barcodes')),what='',quiet=TRUE)
  donors<-fread(get_path('donors'));stopifnot(!anyDuplicated(donors$barcode))
  raw<-as(Matrix::readMM(gzfile(get_path('matrix'))),'CsparseMatrix')
  stopifnot(nrow(raw)==nrow(features),ncol(raw)==length(barcodes))
  keys<-make.unique(gsub('_','-',features[[2]],fixed=TRUE));rownames(raw)<-keys;colnames(raw)<-barcodes
  feature_map<-data.table(feature_key=keys,gene_id=features[[1]],gene_symbol=features[[2]])
  idx<-match(donors$barcode,barcodes);stopifnot(!anyNA(idx))
  toc<-raw[,idx,drop=FALSE];stopifnot(all(Matrix::colSums(toc)>0))
  stage(paste0('preliminary_clustering_for_SoupX ',gsm,' cells=',ncol(toc)))
  tmp<-CreateSeuratObject(toc,project=gsm,min.cells=0,min.features=0)
  tmp<-NormalizeData(tmp,verbose=FALSE);tmp<-FindVariableFeatures(tmp,nfeatures=2000,verbose=FALSE)
  tmp<-ScaleData(tmp,features=VariableFeatures(tmp),verbose=FALSE);tmp<-RunPCA(tmp,npcs=30,seed.use=216187+i,verbose=FALSE)
  tmp<-FindNeighbors(tmp,dims=1:20,verbose=FALSE);tmp<-FindClusters(tmp,resolution=.5,random.seed=216187+i,verbose=FALSE)
  stage(paste0('SoupX_auto_estimation ',gsm))
  sc<-SoupChannel(raw,toc,calcSoupProfile=TRUE);sc<-setClusters(sc,Idents(tmp));sc<-autoEstCont(sc,doPlot=FALSE)
  corrected<-adjustCounts(sc,roundToInt=TRUE)
  rho<-sc$metaData$rho;names(rho)<-rownames(sc$metaData)
  obj<-CreateSeuratObject(corrected,project=gsm,min.cells=0,min.features=0)
  obj[['RNA_uncorrected']]<-CreateAssayObject(counts=toc[,colnames(obj),drop=FALSE])
  di<-match(colnames(obj),donors$barcode);stopifnot(!anyNA(di))
  obj$barcode_original<-colnames(obj);obj$GSM<-gsm;obj$sample_id<-paste0('BM_',group,'_123');obj$group<-group
  obj$GEO_label_before_20260308<-groups$legacy[i]
  obj$souporcell_status<-donors$status[di];obj$souporcell_assignment<-as.character(donors$assignment[di])
  obj$souporcell_log_prob_singleton<-donors$log_prob_singleton[di];obj$souporcell_log_prob_doublet<-donors$log_prob_doublet[di]
  obj$patient_id<-ifelse(obj$souporcell_status=='singlet',paste0(gsm,'_genotype_',obj$souporcell_assignment),NA_character_)
  obj$percent.mt<-PercentageFeatureSet(obj,pattern='^MT-');obj$SoupX_rho<-rho[colnames(obj)]
  q<-as.data.table(obj@meta.data,keep.rownames='barcode');q[,pass_demux:=souporcell_status=='singlet']
  q[,pass_expression_QC:=nFeature_RNA>=200&percent.mt<=15];q[,retained:=pass_demux&pass_expression_QC]
  fwrite(q,paste0('data/processed/',gsm,'_cell_QC.tsv.gz'),sep='\t')
  summaries[[gsm]]<-data.table(GSM=gsm,group=group,raw_matrix_columns=ncol(raw),genotyped_cell_records=nrow(donors),demux_singlets=sum(q$pass_demux),demux_doublets=sum(q$souporcell_status=='doublet'),demux_unassigned=sum(q$souporcell_status=='unassigned'),retained=sum(q$retained),SoupX_rho=median(rho))
  donor_tables[[gsm]]<-as.data.frame(donors);soup_profiles[[gsm]]<-list(rho=median(rho),soup_profile=sc$soupProfile)
  obj<-subset(obj,cells=q[retained==TRUE,barcode]);objects[[gsm]]<-obj
  assign('GSE216187_QC_samples',objects,envir=.GlobalEnv)
  rm(raw,toc,corrected,tmp,sc,obj,barcodes);gc()
  stage(paste0('pool_complete ',gsm,' retained=',sum(q$retained)))
 }
 stage('merge_corrected_RNA')
 obj<-merge(objects[[1]],y=objects[-1],add.cell.ids=names(objects),merge.data=FALSE,project='GSE216187')
 obj$group<-factor(obj$group,levels=c('CR','REL'));obj$patient_id<-factor(obj$patient_id)
 DefaultAssay(obj)<-'RNA';obj<-NormalizeData(obj,verbose=FALSE)
 stage('SCTransform_regress_mito')
 obj<-SCTransform(obj,assay='RNA',new.assay.name='SCT',vars.to.regress='percent.mt',variable.features.n=3000,return.only.var.genes=TRUE,seed.use=216187,verbose=TRUE)
 gene_ids<-sub('\\..*$','',feature_map$gene_id)
 available<-intersect(gene_ids,AnnotationDbi::keys(org.Hs.eg.db,keytype='ENSEMBL'))
 chr<-AnnotationDbi::select(org.Hs.eg.db,keys=available,keytype='ENSEMBL',columns='CHR')
 sex_ids<-unique(chr$ENSEMBL[chr$CHR%in%c('X','Y')]);sex_features<-feature_map$feature_key[gene_ids%in%sex_ids]
 vf<-VariableFeatures(obj);technical<-grepl('^(RP[SL]|MT-)',vf)|vf%in%sex_features
 VariableFeatures(obj)<-vf[!technical];VariableFeatures(obj[['RNA']])<-vf[!technical]
 stage(paste0('PCA_50_components features=',length(VariableFeatures(obj))))
 obj<-RunPCA(obj,assay='SCT',features=VariableFeatures(obj),npcs=50,seed.use=216187,verbose=FALSE)
 choice<-read_json('resources/analysis_choice.json',simplifyVector=TRUE)
 reduction<-'pca'
 if(isTRUE(choice$harmony_by_patient)){
  stage('Harmony_by_souporcell_patient_id')
  obj<-harmony::RunHarmony(obj,group.by.vars='patient_id',reduction='pca',dims.use=1:50,max.iter.harmony=30,project.dim=FALSE)
  reduction<-'harmony'
 }
 stage('neighbors_Louvain_clustering')
 obj<-FindNeighbors(obj,reduction=reduction,dims=1:50,graph.name=c('analysis_nn','analysis_snn'),verbose=FALSE)
 obj<-FindClusters(obj,graph.name='analysis_snn',resolution=c(.5,1,1.5),random.seed=216187,verbose=FALSE)
 obj$seurat_clusters<-obj$analysis_snn_res.1.5;Idents(obj)<-'seurat_clusters'
 stage('UMAP_50_components')
 obj<-RunUMAP(obj,reduction=reduction,dims=1:50,n.neighbors=30,min.dist=.3,seed.use=216187,return.model=TRUE,verbose=TRUE)
 DefaultAssay(obj)<-'RNA'
 qc<-rbindlist(summaries);fwrite(qc,'data/processed/QC_summary.tsv',sep='\t')
 obj@misc$processing<-list(dataset='GSE216187',GEO_identity_correction=list(date='2026-03-08',GSM6660161=list(current='REL',before='CR'),GSM6660162=list(current='CR',before='REL'),file_contents='unchanged'),demultiplexing='Author-provided souporcell singlets; patient identifiers are pool-prefixed genotype IDs',QC=list(min_genes=200,max_mito_pct=15),ambient_RNA=list(method='SoupX',version=as.character(packageVersion('SoupX')),profiles=soup_profiles),normalization='SoupX-corrected RNA counts; RNA LogNormalize for expression; SCT regresses percent.mt',HVG=list(candidate_count=3000,excluded=vf[technical],used=vf[!technical],excluded_classes='ribosomal, mitochondrial, X/Y chromosome genes'),chromosome_annotation=as.character(packageVersion('org.Hs.eg.db')),PCA_components=50,embedding_reduction=reduction,cluster_resolutions=c(.5,1,1.5),default_resolution=1.5,seed=216187,original_donor_tables=donor_tables,feature_map=as.data.frame(feature_map),QC_summary=as.data.frame(qc))
 export_profiles(obj)
 u<-Embeddings(obj,'umap');plot<-data.table(cell_id=colnames(obj),UMAP_1=u[,1],UMAP_2=u[,2],GSM=obj$GSM,group=as.character(obj$group),patient_id=as.character(obj$patient_id),cluster=as.character(obj$seurat_clusters))
 fwrite(plot,'data/plot_ready/UMAP_cells.tsv.gz',sep='\t');fwrite(as.data.table(obj@meta.data,keep.rownames='cell_id'),'data/processed/cell_metadata.tsv.gz',sep='\t')
 write_json(list(status='analysis_ready',cells=ncol(obj),RNA_genes=nrow(obj[['RNA']]),assays=Assays(obj),reductions=Reductions(obj),clusters=nlevels(Idents(obj)),patients=nlevels(obj$patient_id),group_counts=as.list(table(obj$group)),RNA_group_labels='Current GEO corrected labels, 2026-03-08',embedding_reduction=reduction,time=format(Sys.time(),'%FT%T%z')),'logs/analysis_ready.json',auto_unbox=TRUE,pretty=TRUE)
 stage('SC_RNA_ANALYSIS_READY_IN_MEMORY')
 obj
}

stage<-function(message){cat(format(Sys.time(),'%FT%T%z'),message,'\n');flush.console();cat(toJSON(list(time=format(Sys.time(),'%FT%T%z'),stage=message),auto_unbox=TRUE),'\n',file='logs/analysis_progress.jsonl',append=TRUE)}
export_profiles<-function(obj){
 markers<-c('CD34','GATA2','HLF','MEIS1','SPINK2','KIT','MPO','CTSG','AZU1','ELANE','LYZ','ITGAM','LST1','S100A8','IRF7','TCF4','CLEC4C','FCER1A','CD1C','DNTT','VPREB1','CD79A','MS4A1','PPBP','PF4','HBB','GATA1','KLF1','CD3D','CD3E','CD3G','TRAC','TRBC1','TRBC2','CD4','CD8A','CD8B','CCR7','SELL','TCF7','LEF1','IL7R','LTB','ZNF683','GZMK','CD160','NKG7','PRF1','GZMB','GZMH','GNLY','FGFBP2','ADGRG1','TRDC','TRGC1','KLRB1','SLC4A10','NCR3','FOXP3','IL2RA','CTLA4','RORC','CCR6','ISG15','IFI6','IFIT1','MKI67','TOP2A')
 groups<-factor(as.character(obj$seurat_clusters));n<-tabulate(groups)
 design<-sparseMatrix(i=seq_along(groups),j=as.integer(groups),x=1,dims=c(length(groups),nlevels(groups)))
 m<-GetAssayData(obj,assay='RNA',slot='data')[intersect(markers,rownames(obj)),,drop=FALSE]
 avg<-sweep(as.matrix(m%*%design),2,n,'/');pct<-sweep(as.matrix((m>0)%*%design),2,n,'/')*100
 z<-data.table(gene=rep(rownames(m),nlevels(groups)),cluster=rep(levels(groups),each=nrow(m)),mean_log_expression=as.vector(avg),pct_expressing=as.vector(pct))
 fwrite(z,'data/processed/cluster_marker_profiles.tsv',sep='\t')
}
