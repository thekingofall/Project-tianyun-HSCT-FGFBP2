suppressPackageStartupMessages({library(scRepertoire);library(digest)})
source("code/15_memory_bounded_verification.R")

integrate_gse222633_TCR<-function(obj){
 plan<-fread('data/processed/sample_plan.tsv');meta<-as.data.table(obj@meta.data,keep.rownames='cell_id')
 combined<-list();raw_contigs<-list();coverage<-list();primary_all<-list();all_chain_counts<-list()
 baseline<-list(cells=colnames(obj),genes=rownames(obj),RNA_counts=hash_matrix_blocks(GetAssayData(obj,slot='counts')),umap=digest(Embeddings(obj,'umap'),algo='sha256'))
 for(i in seq_len(nrow(plan))){
  s<-plan$sample[i];stage(paste0('TCR_process ',s))
  c<-fread(plan$TCR_file[i],na.strings=c('','None','NA'));c[,source_sample:=s]
  for(z in c('is_cell','high_confidence','full_length','productive'))c[,(z):=bool(get(z))]
  stopifnot(all(c('barcode','chain','cdr3','cdr3_nt','umis','reads','contig_id')%in%names(c)))
  c[,cell_id:=paste0(s,'_',barcode)];raw_contigs[[s]]<-as.data.frame(c)
  c[,valid:=is_cell&high_confidence&full_length&productive&chain%in%c('TRA','TRB','TRG','TRD')&!is.na(cdr3)&grepl('^[ACDEFGHIKLMNPQRSTVWY]{5,}$',cdr3)&!is.na(cdr3_nt)&nzchar(cdr3_nt)]
  valid<-c[valid==TRUE&cell_id%in%colnames(obj)]
  setorder(valid,barcode,chain,-umis,-reads,contig_id)
  unique_chains<-unique(valid,by=c('barcode','chain','cdr3_nt'))
  primary<-unique_chains[,head(.SD,1),by=.(barcode,chain)]
  setcolorder(primary,names(valid))
  primary_all[[s]]<-primary
  chain_counts<-unique_chains[,.(n_chains=.N,all_cdr3_aa=paste(sort(unique(cdr3)),collapse=';'),all_cdr3_nt=paste(sort(unique(cdr3_nt)),collapse=';')),by=.(cell_id,chain)]
  all_chain_counts[[s]]<-chain_counts
  ab<-primary[chain%in%c('TRA','TRB')]
  if(nrow(ab)){
   co<-combineTCR(list(as.data.frame(ab)),samples=s,cells='T-AB',removeNA=FALSE,removeMulti=FALSE,filterMulti=FALSE)
   if(length(co))combined[[s]]<-co[[1]]
  }
  cells<-meta[sample==s,cell_id]
  coverage[[s]]<-data.table(sample=s,RNA_QC_cells=length(cells),TCR_source_barcodes=uniqueN(c$barcode),source_TCR_RNA_overlap=uniqueN(c[cell_id%in%cells,barcode]),valid_any_chain=uniqueN(valid$barcode))
 }
 primary<-rbindlist(primary_all,fill=TRUE);chain_counts<-rbindlist(all_chain_counts,fill=TRUE)
 extra<-data.table(cell_id=colnames(obj))
 for(ch in c('TRA','TRB','TRG','TRD')){
  cc<-chain_counts[chain==ch];x<-primary[chain==ch];idx<-match(extra$cell_id,x$cell_id);ci<-match(extra$cell_id,cc$cell_id)
  extra[,(paste0('TCR_n_',ch)):=ifelse(is.na(ci),0L,cc$n_chains[ci])]
  extra[,(paste0('TCR_valid_',ch)):=!is.na(idx)]
  for(field in c('cdr3','cdr3_nt','v_gene','d_gene','j_gene','c_gene','umis','reads')){
   key<-switch(field,cdr3=paste0('TCR_cdr3_aa_',ch),cdr3_nt=paste0('TCR_cdr3_nt_',ch),paste0('TCR_',field,'_',ch))
   extra[,(key):=x[[field]][idx]]
  }
  extra[,(paste0('TCR_all_cdr3_aa_',ch)):=cc$all_cdr3_aa[ci]]
  extra[,(paste0('TCR_all_cdr3_nt_',ch)):=cc$all_cdr3_nt[ci]]
 }
 extra[,`:=`(TCR_paired_TRA_TRB=TCR_valid_TRA&TCR_valid_TRB,TCR_paired_TRG_TRD=TCR_valid_TRG&TCR_valid_TRD,TCR_any_valid=TCR_valid_TRA|TCR_valid_TRB|TCR_valid_TRG|TCR_valid_TRD,TCR_multichain=TCR_n_TRA>1|TCR_n_TRB>1|TCR_n_TRG>1|TCR_n_TRD>1)]
 extra[,TCR:=ifelse(TCR_paired_TRA_TRB|TCR_paired_TRG_TRD,'TCR+','TCR-')]
 extra[,TCR_any:=ifelse(TCR_any_valid,'TCR+','TCR-')]
 extra[,TCR_unambiguous_TRA_TRB:=TCR_n_TRA==1&TCR_n_TRB==1]
 pair_id<-as.character(obj$pair_id)
 extra[,TCR_pair_clone_id_aa:=ifelse(TCR_paired_TRA_TRB,paste(pair_id,TCR_cdr3_aa_TRA,TCR_cdr3_aa_TRB,sep='::'),NA_character_)]
 extra[,TCR_pair_clone_id_nt:=ifelse(TCR_paired_TRA_TRB,paste(pair_id,TCR_cdr3_nt_TRA,TCR_cdr3_nt_TRB,sep='::'),NA_character_)]
 if(any(extra$TCR_paired_TRG_TRD))stage(paste0('Gamma-delta paired annotations retained: ',sum(extra$TCR_paired_TRG_TRD),'; native clone tables use the study alpha-beta libraries'))
 native_cols<-c('barcode','CTgene','CTnt','CTaa','CTstrict','Frequency','cloneType')
 stopifnot(!any(native_cols%in%names(obj@meta.data)),all(unlist(lapply(combined,function(x)x$barcode))%in%colnames(obj)))
 obj<-combineExpression(combined,obj,cloneCall='aa',chain='both',group.by='none',proportion=TRUE,filterNA=FALSE)
 dd<-as.data.frame(extra[,setdiff(names(extra),'cell_id'),with=FALSE]);rownames(dd)<-extra$cell_id;obj<-AddMetaData(obj,dd)
 clone_cells<-data.table(cell_id=colnames(obj),sample=obj$sample,pair_id=as.character(obj$pair_id),role=as.character(obj$role),timepoint=as.character(obj$timepoint),CTaa=obj$CTaa,pair_clone_id=obj$TCR_pair_clone_id_aa)
 clone_cells[!is.na(CTaa),clonalFrequency:=.N,by=.(sample,CTaa)]
 clone_cells[!is.na(CTaa),clonalProportion:=clonalFrequency/.N,by=sample]
 obj$clonalFrequency<-clone_cells$clonalFrequency;obj$clonalProportion<-clone_cells$clonalProportion
 obj$cloneSize<-cut(obj$clonalFrequency,breaks=c(0,1,5,20,100,Inf),labels=c('Singleton (1)','Small (2-5)','Medium (6-20)','Large (21-100)','Hyperexpanded (>100)'))
 paired_cells<-clone_cells[!is.na(pair_clone_id)]
 pair_table<-paired_cells[,.(n_cells=.N,seen_in_donor=any(role=='Donor'),seen_in_recipient=any(role=='Recipient'),n_samples=uniqueN(sample),timepoints=paste(sort(unique(timepoint)),collapse=';')),by=.(pair_id,pair_clone_id)]
 pair_table[,persisting:=seen_in_donor&seen_in_recipient]
 obj$TCR_persisting<-obj$TCR_pair_clone_id_aa%in%pair_table[persisting==TRUE,pair_clone_id]
 coverage<-rbindlist(coverage)
 coverage<-merge(coverage,extra[,.(cell_id,TCR_paired_TRA_TRB,TCR_multichain)][,sample:=obj$sample][,.(paired_TRA_TRB=sum(TCR_paired_TRA_TRB),multi_chain_cells=sum(TCR_multichain)),by=sample],by='sample',all.x=TRUE)
 stopifnot(identical(baseline$cells,colnames(obj)),identical(baseline$genes,rownames(obj)),identical(baseline$RNA_counts,hash_matrix_blocks(GetAssayData(obj,slot='counts'))),identical(baseline$umap,digest(Embeddings(obj,'umap'),algo='sha256')))
 stopifnot(all(abs(obj$Frequency-obj$clonalProportion)<1e-10,na.rm=TRUE))
 obj@misc$TCR<-list(combined=combined,all_source_contigs=raw_contigs,selected_primary_contigs=as.data.frame(primary),coverage=as.data.frame(coverage),paired_clone_summary=as.data.frame(pair_table),criteria='is_cell & high_confidence & full_length & productive; complete valid CDR3 amino acid and nucleotide sequences; exact sample+barcode',primary_chain_selection='Per barcode and chain: highest UMI, then reads, then contig_id; all secondary sequences and raw contigs retained',TCR_definition='At least one valid TRA and one valid TRB in same cell; paired TRG/TRD flags are also retained',clone_definition='Native CTaa for general repertoire; pair_id+paired primary CDR3 AA/NT for donor-recipient tracking',persistence_definition='Paired AA clonotype seen in at least one donor and at least one recipient sample of same pair',scRepertoire_version=as.character(packageVersion('scRepertoire')))
 fwrite(coverage,'data/processed/TCR_coverage.tsv',sep='\t');fwrite(pair_table,'data/processed/paired_clone_tracking.tsv.gz',sep='\t')
 fwrite(as.data.table(obj@meta.data,keep.rownames='cell_id'),'data/processed/cell_metadata.tsv.gz',sep='\t')
 u<-Embeddings(obj,'umap');plot<-data.table(cell_id=colnames(obj),UMAP_1=u[,1],UMAP_2=u[,2],sample=obj$sample,pair_id=as.character(obj$pair_id),role=as.character(obj$role),timepoint=as.character(obj$timepoint),cluster=as.character(obj$seurat_clusters),TCR=obj$TCR,TCR_any=obj$TCR_any,TCR_paired_TRA_TRB=obj$TCR_paired_TRA_TRB,TCR_multichain=obj$TCR_multichain,TCR_persisting=obj$TCR_persisting)
 fwrite(plot,'data/plot_ready/RNA_TCR_UMAP_cells.tsv.gz',sep='\t')
 summary<-list(status='analysis_ready',cells=ncol(obj),genes=nrow(obj),assays=Assays(obj),reductions=Reductions(obj),clusters=nlevels(Idents(obj)),samples=nrow(plan),pairs=nlevels(obj$pair_id),TCR_any=sum(obj$TCR_any_valid),TCR_paired=sum(obj$TCR=='TCR+'),TCR_unambiguous_paired=sum(obj$TCR_unambiguous_TRA_TRB),TCR_multichain=sum(obj$TCR_multichain),persisting_paired_clonotypes=sum(pair_table$persisting),persisting_cells=sum(obj$TCR_persisting),time=format(Sys.time(),'%FT%T%z'))
 write_json(summary,'logs/analysis_ready.json',auto_unbox=TRUE,pretty=TRUE);print(summary)
 stage('RNA_TCR_ANALYSIS_READY_IN_MEMORY')
 obj
}
bool<-function(x)tolower(as.character(x))%in%c('true','t','1')
