# Run in the existing R session; source object and original files are retained.
suppressPackageStartupMessages({library(Seurat);library(Matrix);library(data.table);library(jsonlite)})
local({
root<-'/data5/tiany/GSE222633/GSE222633_Teff_longitudinal_v20261031_01';setwd(root)
con<-file('logs/02_extract.log','wt');sink(con,split=TRUE);sink(con,type='message');on.exit({sink(type='message');sink();close(con)},add=TRUE)
input<-'/data5/tiany/GSE222633/CD8.rds';before<-file.info(input)[,c('size','mtime')]
if(!exists('CD8_candidate_GSE222633',envir=.GlobalEnv))assign('CD8_candidate_GSE222633',readRDS(input),envir=.GlobalEnv)
ob<-get('CD8_candidate_GSE222633',envir=.GlobalEnv)
md<-as.data.table(ob@meta.data,keep.rownames='cell');stopifnot(nrow(md)==ncol(ob),identical(md$cell,colnames(ob)),all(c('Celltype','Pair','sample','timepoint','role')%in%names(md)))
stopifnot('CD8 Teff'%in%md$Celltype)
md[,`:=`(Pair=as.character(Pair),timepoint=as.character(timepoint),Celltype=as.character(Celltype))]
stopifnot(all(md$Pair%in%LETTERS[1:5]),all(md$timepoint%in%c('Pre_GCSF','Post_GCSF','D90','D180')))
counts<-GetAssayData(ob,assay='RNA',slot='counts');den<-Matrix::colSums(counts);stopifnot(all(den>0))
gsets<-fread('resources/state_gene_sets.tsv');genes<-unique(c('FGFBP2','ADGRG1',gsets$gene));missing<-setdiff(genes,rownames(counts));stopifnot(all(c('FGFBP2','ADGRG1')%in%rownames(counts)))
fwrite(data.table(gene=missing),'reports/missing_genes.csv')
genes<-intersect(genes,rownames(counts));norm<-counts[genes,,drop=FALSE]%*%Diagonal(x=10000/den);rownames(norm)<-genes;colnames(norm)<-colnames(counts);lg<-log1p(norm)
md[,`:=`(FGFBP2=as.numeric(lg['FGFBP2',]),ADGRG1=as.numeric(lg['ADGRG1',]),FGFBP2_positive=as.numeric(counts['FGFBP2',])>0,ADGRG1_positive=as.numeric(counts['ADGRG1',])>0)]
md[,double_positive:=FGFBP2_positive & ADGRG1_positive]
for(pr in unique(gsets$program)){gs<-intersect(gsets[program==pr,gene],genes);stopifnot(length(gs)>=4,!any(c('FGFBP2','ADGRG1')%in%gs));md[,(pr):=as.numeric(Matrix::colMeans(lg[gs,,drop=FALSE]))]}
sm<-md[,.(n_CD8=.N,n_Teff=sum(Celltype=='CD8 Teff')),by=.(sample,Pair,participant_id,role,timepoint)];sm[,Teff_percent:=100*n_Teff/n_CD8];fwrite(sm,'data/plot_ready/sample_composition.csv')
te<-md[Celltype=='CD8 Teff'];cat('CD8',nrow(md),'Teff',nrow(te),'\n')
ex<-rbindlist(lapply(c('FGFBP2','ADGRG1'),function(g){te[,.(gene=g,n_cells=.N,mean_log_expression=mean(get(g)),log1p_mean_normalized=log1p(mean(expm1(get(g)))),positive_percent=100*mean(get(paste0(g,'_positive')))),by=.(sample,Pair,participant_id,role,timepoint)]}))
fwrite(ex,'data/plot_ready/Teff_gene_expression.csv')
pos<-te[,.(n_cells=.N,FGFBP2=100*mean(FGFBP2_positive),ADGRG1=100*mean(ADGRG1_positive),Double_positive=100*mean(double_positive)),by=.(sample,Pair,participant_id,role,timepoint)];fwrite(pos,'data/plot_ready/Teff_positive_percent.csv')
sc<-te[,lapply(.SD,mean),by=.(sample,Pair,participant_id,role,timepoint),.SDcols=unique(gsets$program)];fwrite(sc,'data/plot_ready/Teff_functional_scores.csv')
# Strict paired identity is scoped to each biological donor-recipient pair.
chaincols<-c('TCR_cdr3_nt_TRA','TCR_cdr3_nt_TRB','TCR_v_gene_TRA','TCR_j_gene_TRA','TCR_v_gene_TRB','TCR_j_gene_TRB')
stopifnot(all(c(chaincols,'TCR_unambiguous_TRA_TRB','TCR_valid_TRA','TCR_valid_TRB')%in%names(md)))
md[,eligible:=TCR_unambiguous_TRA_TRB%in%TRUE & TCR_valid_TRA%in%TRUE & TCR_valid_TRB%in%TRUE]
for(v in chaincols)md[,eligible:=eligible & !is.na(get(v)) & nzchar(trimws(as.character(get(v)))) & !get(v)%in%c('NA','None','none')]
md[,clone_key:=NA_character_];md[eligible==TRUE,clone_key:=do.call(paste,c(.SD,sep='::')),.SDcols=c('Pair',chaincols)]
rec<-unique(md[eligible==TRUE & role=='Recipient',clone_key]);md[,recipient_detected:=eligible & clone_key%in%rec]
dt<-md[Celltype=='CD8 Teff' & role=='Donor' & eligible==TRUE];dt[,Expression_group:=ifelse(double_positive,'Double positive','Other Teff')]
ps<-dt[,.(eligible_cells=.N,recipient_detected_cells=sum(recipient_detected),recipient_detected_percent=100*mean(recipient_detected),clonotypes=uniqueN(clone_key)),by=.(sample,Pair,timepoint,Expression_group)];fwrite(ps,'data/plot_ready/donor_Teff_recipient_detection.csv')
fwrite(md[,.(n_CD8=.N,paired_cells=sum(eligible),n_Teff=sum(Celltype=='CD8 Teff'),Teff_paired_cells=sum(eligible & Celltype=='CD8 Teff')),by=.(sample,Pair,timepoint)],'reports/TCR_coverage.csv')
fwrite(te[,.(cell,sample,Pair,participant_id,role,timepoint,FGFBP2,ADGRG1,FGFBP2_positive,ADGRG1_positive,double_positive,Cytotoxicity,Exhaustion,Naive_memory)],'data/plot_ready/Teff_cell_values.tsv',sep='\t')
rank90<-ex[timepoint=='D90',.(Pair,gene,n_cells,mean_log_expression,positive_percent)];rank90[,rank_low_to_high:=frank(mean_log_expression,ties.method='min'),by=gene];fwrite(rank90,'reports/D90_expression_ranks.csv')
dd<-dcast(ex[Pair=='D'],gene~timepoint,value.var='mean_log_expression');if(all(c('D90','D180')%in%names(dd)))dd[,change_D180_minus_D90:=D180-D90];fwrite(dd,'reports/Pair_D_expression_change.csv')
stopifnot(all(ex$n_cells>0),all(ex$positive_percent>=0 & ex$positive_percent<=100),identical(before,file.info(input)[,c('size','mtime')]))
write_json(list(CD8_cells=nrow(md),Teff_cells=nrow(te),n_samples=nrow(sm),n_pairs=uniqueN(sm$Pair),source_unchanged=TRUE,assay='RNA counts',normalization='log1p(counts / total RNA UMI * 10000)',primary_summary='mean per-cell log1p normalized expression',positive='raw count > 0',source_mtime=as.character(before$mtime)), 'reports/extraction_QA.json',pretty=TRUE,auto_unbox=TRUE)
writeLines(capture.output(sessionInfo()),'logs/sessionInfo.txt');print(ex[Pair=='D']);print(rank90);cat('DONE_EXTRACT\n')
source('code/03_plot.R')
})
