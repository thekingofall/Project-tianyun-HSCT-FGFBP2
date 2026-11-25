suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(DoubletFinder);library(harmony);library(jsonlite)})
setDTthreads(2);set.seed(222633);options(warn=1,future.globals.maxSize=16*1024^3)
future::plan('sequential')

process_gse222633_RNA<-function(){
 plan<-fread('data/processed/sample_plan.tsv')
 samples<-if(exists('GSE222633_QC_samples',envir=.GlobalEnv))get('GSE222633_QC_samples',envir=.GlobalEnv)else list()
 qc<-if(exists('GSE222633_QC_audit',envir=.GlobalEnv))get('GSE222633_QC_audit',envir=.GlobalEnv)else list()
 pending<-which(!plan$sample%in%names(samples));last_wait<-0
 dir.create('data/processed/QC',showWarnings=FALSE)
 while(length(pending)){
  ready<-pending[vapply(pending,function(i)file_ready(plan$RNA_file[i])&&file_ready(plan$TCR_file[i]),TRUE)]
  if(!length(ready)){
   if(as.numeric(Sys.time())-last_wait>60){stage(paste0('waiting_for_downloads: ',length(pending),' samples remaining'));last_wait<-as.numeric(Sys.time())}
   Sys.sleep(15);next
  }
  for(i in ready){
   set.seed(222633+i);s<-plan$sample[i];stage(paste0('read_and_QC ',s))
   raw<-Read10X_h5(plan$RNA_file[i],use.names=TRUE,unique.features=TRUE)
   stopifnot(is.list(raw),'Gene Expression'%in%names(raw))
   obj<-CreateSeuratObject(raw[['Gene Expression']],project=s,min.cells=0,min.features=0)
   obj$percent.mt<-PercentageFeatureSet(obj,pattern='^MT-')
   obj$sample<-s;obj$pair_id<-as.character(plan$pair_id[i]);obj$participant_id<-plan$participant_id[i];obj$role<-plan$role[i];obj$timepoint<-plan$timepoint[i]
   qcmd<-as.data.table(obj@meta.data,keep.rownames='barcode');qcmd[,gene_mt_pass:=nFeature_RNA>=250&nFeature_RNA<=5000&percent.mt<10]
   n_raw<-ncol(obj);obj<-subset(obj,cells=qcmd[gene_mt_pass==TRUE,barcode]);n_gene<-ncol(obj)
   if('Antibody Capture'%in%names(raw)){
    adt<-raw[['Antibody Capture']][,colnames(obj),drop=FALSE];ht<-grepl('Hashtag[0-9]',rownames(adt),ignore.case=TRUE)
    if(any(!ht)){obj[['ADT']]<-CreateAssayObject(counts=adt[!ht,,drop=FALSE]);obj<-NormalizeData(obj,assay='ADT',normalization.method='CLR',margin=2,verbose=FALSE)}
    if(any(ht)){obj[['HTO']]<-CreateAssayObject(counts=adt[ht,,drop=FALSE]);obj<-NormalizeData(obj,assay='HTO',normalization.method='CLR',margin=2,verbose=FALSE)}
   }
   rm(raw);gc()
   obj<-NormalizeData(obj,verbose=FALSE);obj<-FindVariableFeatures(obj,nfeatures=2000,verbose=FALSE)
   obj<-ScaleData(obj,features=VariableFeatures(obj),verbose=FALSE);obj<-RunPCA(obj,npcs=30,verbose=FALSE,seed.use=222633+i)
   obj<-FindNeighbors(obj,dims=1:20,verbose=FALSE);obj<-FindClusters(obj,resolution=.5,verbose=FALSE,random.seed=222633+i)
   homotypic<-DoubletFinder::modelHomotypic(as.character(obj$seurat_clusters));n_exp<-round(.05*ncol(obj));n_adj<-round(n_exp*(1-homotypic))
   stage(paste0('DoubletFinder ',s,' n=',ncol(obj),' expected_adjusted=',n_adj))
   fun<-get0('doubletFinder_v3',envir=asNamespace('DoubletFinder'))
   if(is.null(fun))fun<-get('doubletFinder',envir=asNamespace('DoubletFinder'))
   result<-fun(obj,PCs=1:20,pN=.25,pK=.09,nExp=n_adj,reuse.pANN=FALSE,sct=FALSE)
   fc<-tail(grep('^DF.classifications',names(result@meta.data),value=TRUE),1)
   pc<-tail(grep('^pANN',names(result@meta.data),value=TRUE),1)
   obj$doublet_class<-result@meta.data[[fc]];obj$doublet_pANN<-result@meta.data[[pc]]
   audit<-data.table(barcode=colnames(obj),doublet_class=obj$doublet_class,doublet_pANN=obj$doublet_pANN)
   full_audit<-merge(qcmd,audit,by='barcode',all.x=TRUE,sort=FALSE)
   fwrite(full_audit,paste0('data/processed/QC/',s,'_cell_QC.tsv.gz'),sep='\t')
   obj<-subset(obj,subset=doublet_class=='Singlet')
   qc[[s]]<-data.table(sample=s,pair_id=plan$pair_id[i],role=plan$role[i],timepoint=plan$timepoint[i],raw_cells=n_raw,gene_mt_pass=n_gene,predicted_doublets=n_gene-ncol(obj),retained_cells=ncol(obj),homotypic_fraction=homotypic,expected_doublets=n_exp,adjusted_expected_doublets=n_adj)
   obj$QC_precluster<-as.character(obj$seurat_clusters);obj$RNA_snn_res.0.5<-NULL;obj$seurat_clusters<-NULL
   obj<-DietSeurat(obj,counts=TRUE,data=TRUE,scale.data=FALSE,assays=Assays(obj),dimreducs=NULL,graphs=NULL)
   samples[[s]]<-obj
   assign('GSE222633_QC_samples',samples,envir=.GlobalEnv);assign('GSE222633_QC_audit',qc,envir=.GlobalEnv)
   rm(result,obj,qcmd,audit,full_audit);gc()
   pending<-setdiff(pending,i);stage(paste0('sample_complete ',s,' remaining=',length(pending)))
  }
 }
 stage('merge_QC_samples')
 names_order<-plan$sample;samples<-samples[names_order]
 obj<-merge(samples[[1]],y=samples[-1],add.cell.ids=names(samples),merge.data=TRUE,project='GSE222633')
 rm(samples);gc()
 obj$pair_id<-factor(obj$pair_id);obj$role<-factor(obj$role,levels=c('Donor','Recipient'));obj$timepoint<-factor(obj$timepoint,levels=c('Pre_GCSF','Post_GCSF','D90','D180'))
 DefaultAssay(obj)<-'RNA'
 stage('global_normalization_and_HVG')
 obj<-NormalizeData(obj,verbose=FALSE);obj<-FindVariableFeatures(obj,nfeatures=3000,verbose=FALSE)
 hv<-VariableFeatures(obj);exclude_vdj<-grepl('^(TR[ABDG][VDJ]|IG[HKL][VJ]|IGHD[0-9])',hv)
 VariableFeatures(obj)<-hv[!exclude_vdj]
 obj<-ScaleData(obj,features=VariableFeatures(obj),verbose=FALSE)
 stage('PCA_50_components');obj<-RunPCA(obj,npcs=50,verbose=FALSE,seed.use=222633)
 stage('Harmony_by_donor_recipient_pair')
 obj<-harmony::RunHarmony(obj,group.by.vars='pair_id',reduction='pca',dims.use=1:30,theta=2,max.iter.harmony=30,project.dim=FALSE,verbose=TRUE)
 stage('neighbors_and_clustering')
 obj<-FindNeighbors(obj,reduction='harmony',dims=1:30,k.param=20,graph.name=c('harmony_nn','harmony_snn'),verbose=FALSE)
 obj<-FindClusters(obj,graph.name='harmony_snn',resolution=c(.3,.6,1),random.seed=222633,verbose=FALSE)
 obj$seurat_clusters<-obj$harmony_snn_res.0.6;Idents(obj)<-'seurat_clusters'
 stage('UMAP_from_Harmony')
 obj<-RunUMAP(obj,reduction='harmony',dims=1:30,n.neighbors=30,min.dist=.3,seed.use=222633,return.model=TRUE,verbose=TRUE)
 qc_table<-rbindlist(qc);fwrite(qc_table,'data/processed/QC_summary.tsv',sep='\t')
 obj@misc$RNA_processing<-list(dataset='GSE222633',QC=list(min_genes=250,max_genes=5000,max_mito_pct=10,doublet_method='DoubletFinder',doublet_rate=.05,doublet_pN=.25,doublet_pK=.09,homotypic_adjustment='Preliminary per-sample RNA clusters'),normalization='LogNormalize 10000',ADT='CLR margin 2; hashtag features retained separately in HTO',hvg_requested=3000,hvg_used=length(VariableFeatures(obj)),excluded_variable_VDJ_genes=hv[exclude_vdj],PCA_components=50,Harmony=list(group='pair_id',theta=2,dims=1:30,max_iter=30),neighbor_k=20,cluster_resolutions=c(.3,.6,1),active_resolution=.6,UMAP=list(reduction='harmony',dims=1:30,min_dist=.3,n_neighbors=30),seed=222633,reference_prediction_filter='Not applied; Harmony integration selected by user',sample_manifest=as.data.frame(plan),QC_summary=as.data.frame(qc_table))
 export_cluster_profiles(obj)
 stage(paste0('RNA_processing_complete cells=',ncol(obj),' genes=',nrow(obj)))
 obj
}

file_ready<-function(path)file.exists(path)&&!file.exists(paste0(path,'.aria2'))&&unname(file.info(path)$size)>0
stage<-function(message){cat(format(Sys.time(),'%FT%T%z'),message,'\n');flush.console();line<-toJSON(list(time=format(Sys.time(),'%FT%T%z'),stage=message),auto_unbox=TRUE);cat(line,'\n',file='logs/analysis_progress.jsonl',append=TRUE)}
export_cluster_profiles<-function(obj){
 markers<-c('CD3D','CD3E','TRAC','CD4','IL7R','CCR7','LEF1','TCF7','LTB','CD8A','CD8B','GZMK','GZMB','GNLY','NKG7','PRF1','FGFBP2','KLRD1','TRDC','TRGC1','TRGC2','FOXP3','IL2RA','CTLA4','TIGIT','MKI67','TOP2A','STMN1','MS4A1','CD79A','CD79B','CD74','CD14','LYZ','FCGR3A','LST1','S100A8','S100A9','FCER1A','CD1C','CLEC10A','CLEC4C','TCF4','IRF7','MZB1','JCHAIN','XBP1','PPBP','PF4','HBB','HBA1','CD34','GATA2')
 groups<-factor(as.character(obj$seurat_clusters));design<-sparseMatrix(i=seq_along(groups),j=as.integer(groups),x=1,dims=c(length(groups),nlevels(groups)))
 gm<-GetAssayData(obj,slot='data')[intersect(markers,rownames(obj)),,drop=FALSE]
 avg<-as.matrix(gm%*%design);avg<-sweep(avg,2,tabulate(groups),'/');pct<-as.matrix((gm>0)%*%design);pct<-sweep(pct,2,tabulate(groups),'/')*100
 z<-as.data.table(as.table(avg));setnames(z,c('gene','cluster_index','mean_log_expression'));z[,cluster:=levels(groups)[as.integer(cluster_index)]];z[,pct_expressing:=as.vector(pct)]
 fwrite(z,'data/processed/cluster_marker_profiles.tsv',sep='\t')
 if('ADT'%in%Assays(obj)){
  a<-GetAssayData(obj,assay='ADT',slot='data');a<-as.matrix(a%*%design);a<-sweep(a,2,tabulate(groups),'/')
  z<-as.data.table(as.table(a));setnames(z,c('protein','cluster_index','mean_CLR'));z[,cluster:=levels(groups)[as.integer(cluster_index)]];fwrite(z,'data/processed/cluster_ADT_profiles.tsv',sep='\t')
 }
 u<-Embeddings(obj,'umap');d<-data.table(cell_id=colnames(obj),UMAP_1=u[,1],UMAP_2=u[,2],sample=obj$sample,pair_id=as.character(obj$pair_id),role=as.character(obj$role),timepoint=as.character(obj$timepoint),cluster=as.character(obj$seurat_clusters))
 fwrite(d,'data/plot_ready/RNA_UMAP_cells.tsv.gz',sep='\t');fwrite(data.table(PC=seq_along(Stdev(obj,'pca')),stdev=Stdev(obj,'pca')),'data/plot_ready/PCA_stdev.tsv',sep='\t')
}
