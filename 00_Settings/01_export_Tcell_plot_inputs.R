suppressPackageStartupMessages(library(Seurat))
args <- commandArgs(TRUE)
src <- args[1]; out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_table <- function(d,name) {f=gzfile(file.path(out,paste0(name,'.csv.gz')),'wt');write.csv(d,f,row.names=FALSE);close(f)}
genes <- c('CD34','SPINK2','CD3D','CD3E','CD3G','NKG7','GNLY','CD79A','CD19','CD14','VCAN','FCGR3B','S100A8','FCER1A','CLEC4C','IL7R','CCR7','TCF7','SELL','PRF1','GZMH','GZMK','GZMB','FGFBP2','ADGRG1','LAG3','TIGIT','PDCD1','CTLA4','HAVCR2','TOX','TBX21','EOMES','PRDM1','FOXP3','IL2RA','ZNF683','KLRB1','SLC4A10','TRDC','TRGC1','TRGC2','CD8A','CD8B','CD4','MKI67','TOP2A')
x <- readRDS(file.path(src,'Tcell.rds'))
stopifnot(ncol(x)==72884L,length(unique(x$orig.ident))==23L)
md <- x@meta.data;md$cell_id <- rownames(md)
write_table(md,'Tcell_metadata')
u <- as.data.frame(Embeddings(x,'umap'));u$cell_id <- rownames(u);write_table(u,'Tcell_UMAP')
g <- intersect(genes,rownames(x[['RNA']]))
e <- as.data.frame(t(as.matrix(GetAssayData(x,assay='RNA',slot='data')[g,,drop=FALSE])));e$cell_id <- rownames(e);write_table(e,'Tcell_marker_logexpression')
for(v in c('Celltype','Celltype2','condition')) {d=as.data.frame(table(sample=x$orig.ident,label=x[[v]][,1]));write.csv(d,file.path(out,paste0('Tcell_sample_',v,'_counts.csv')),row.names=FALSE)}
rm(x,e);gc()
for(nm in c('CD8cds','CD4cds')) {
 x <- readRDS(file.path(src,paste0(nm,'.rds')))
 md <- Biobase::pData(x);md$cell_id=rownames(md);write_table(md,paste0(nm,'_metadata'))
 coords <- t(x@reducedDimS);d=as.data.frame(coords);d$cell_id=rownames(d);write_table(d,paste0(nm,'_trajectory_coordinates'))
 saveRDS(list(coordinates=x@reducedDimS,branch_points=x@reducedDimK,minSpanningTree=x@minSpanningTree),file.path(out,paste0(nm,'_trajectory_geometry.rds')))
 rm(x);gc()
}
writeLines(capture.output(sessionInfo()),file.path(out,'export_sessionInfo.txt'))
cat('Export finished without refitting or reannotation\n')
