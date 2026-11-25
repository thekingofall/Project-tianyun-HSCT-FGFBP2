library(scRepertoire)
library(Seurat)
library(ggsci)

# TCR导入文件------------------------
GZ0006T <- read.csv("../runs/GZ0006/outs/per_sample_outs/GZ0006/vdj_t/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0006T <- read.csv("../runs/YZ0006/outs/per_sample_outs/YZ0006/vdj_t/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0006aT <- read.csv("../runs/YZ0006a/outs/per_sample_outs/YZ0006a/vdj_t/filtered_contig_annotations.csv",
                     stringsAsFactors = F)
GZ0007T <- read.csv("../runs/GZ0007/outs/per_sample_outs/GZ0007/vdj_t/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0007T <- read.csv("../runs/YZ0007/outs/per_sample_outs/YZ0007/vdj_t/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0007aT <- read.csv("../runs/YZ0007a/outs/per_sample_outs/YZ0007a/vdj_t/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
contig_list <- list(GZ0006T,YZ0006T,YZ0006aT,GZ0007T,YZ0007T,YZ0007aT)
head(contig_list[[1]])
names(contig_list)
#合并数据
combined <- combineTCR(contig_list, samples = c('GZ0006','YZ0006','YZ0006a','GZ0007','YZ0007','YZ0007a'),
                       cells ="T-AB")#

combined <- addVariable(combined, name = "ID", variables = c("Donor","Pre HSCT","Pro HSCT","Donor","Pre HSCT","Pro HSCT"))

head(combined$GZ0007)

#可视化
#1.唯一克隆类型的总数或相对数量
p1 <-quantContig(combined, cloneCall="gene+nt",scale = T)+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))#比例
p2 <- quantContig(combined, cloneCall="gene+nt", scale = F)+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))#个数
p1 + p2# “gene+nt” 使用包含 TCR/Ig + CDR3 区域的核苷酸序列的 VDJC 基因
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = T, 
                                  exportTable = T)#输出
quantContig_output
#2.克隆型的相对分布（丰度）
abundanceContig(combined, cloneCall = "gene", scale = F)+scale_fill_lancet()
#3.CDR3序列的长度分布，"nt"CDR3核苷酸序列"aa"氨基酸序列
lengthContig(combined, cloneCall="aa", chain = "single")+scale_fill_npg()
lengthContig(combined, cloneCall="nt")+scale_fill_npg()
#4.样本之间的克隆类型和动态变化，没有共享则没有连线
compareClonotypes(combined,numbers = 16,samples = c("GZ0006T_Donor", "YZ0006aT_Pro HSCT"),
                  cloneCall="aa", graph = "alluvial") +  scale_fill_npg()
#5.克隆空间稳态
clonalHomeostasis(combined, cloneCall = "gene",
                  cloneTypes = c(Rare = 1e-04,
                                 Small = 0.001,
                                 Medium = 0.01,
                                 Large = 0.1,
                                 Hyperexpanded = 1))
#6.克隆比例，按克隆型的出现频率将其进行排名，1:10表示每个样本中的前10个克隆型
clonalProportion(combined, cloneCall = "nt")
#7. Overlap Analysis 重叠分析
clonalOverlap(combined, cloneCall = "gene+nt",   
              method = "overlap")

clonesizeDistribution(combined, 
                      cloneCall = "gene+nt", 
                      method="ward.D2")
#8.多样性分析
clonalDiversity(combined, 
                cloneCall = "gene+nt",  
                n.boots = 100)

clonalDiversity(combined, 
                cloneCall = "gene+nt",  
                group = "ID",
                n.boots = 100)
#和TCR合并
scRNA <- combineExpression(combined, scRNA, 
                            cloneCall="gene", 
                            # group.by = "Sample", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
head(scRNA1)
DimPlot(scRNA, group.by = "cloneType",label = F)+
  DimPlot(scRNA, group.by = "name",label = T)
#cloneType分布柱形图
seurat_T <- subset(Tcell, tc %in% c("CD4_Naive","CD8_Naive","CD4_CM","CD8_CM","CD4_EM","CD8_EM",
                                               "CD4_TE","CD8_TE","Th","Treg"))
occupiedscRepertoire(seurat_T, x.axis = "tc")
Idents(Tcell) <- 'tc'
combined2 <- expression2List(Tcell,split.by = 'tc')
length(combined2) 
clonalHomeostasis(combined2, 
                  cloneCall = "nt")
clonalDiversity(combined2, 
                cloneCall = "nt")
clonalOverlap(combined2, 
              cloneCall="aa", 
              method="overlap")
# BCR------------------------
GZ0006B <- read.csv("../runs/GZ0006/outs/per_sample_outs/GZ0006/vdj_b/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0006B <- read.csv("../runs/YZ0006/outs/per_sample_outs/YZ0006/vdj_b/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0006aB <- read.csv("../runs/YZ0006a/outs/per_sample_outs/YZ0006a/vdj_b/filtered_contig_annotations.csv",
                     stringsAsFactors = F)
GZ0007B <- read.csv("../runs/GZ0007/outs/per_sample_outs/GZ0007/vdj_b/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0007B <- read.csv("../runs/YZ0007/outs/per_sample_outs/YZ0007/vdj_b/filtered_contig_annotations.csv",
                    stringsAsFactors = F)
YZ0007aB <- read.csv("../runs/YZ0007a/outs/per_sample_outs/YZ0007a/vdj_b/filtered_contig_annotations.csv",
                     stringsAsFactors = F)
contig_list <- list(GZ0006B,YZ0006B,YZ0006aB,GZ0007B,YZ0007B,YZ0007aB)
head(contig_list[[1]])
names(contig_list)
#合并数据
combined <- combineBCR(contig_list, samples = c('GZ0006','YZ0006','YZ0006a','GZ0007','YZ0007','YZ0007a'))#

combined <- addVariable(combined, name = "ID", variables = c("Donor","Pre HSCT","Pro HSCT","Donor","Pre HSCT","Pro HSCT"))

head(combined$GZ0007)
#可视化
#1.唯一克隆类型的总数或相对数量
p1 <-quantContig(combined, cloneCall="gene+nt",scale = T)+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))#比例
p2 <- quantContig(combined, cloneCall="gene+nt", scale = F)+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.85,hjust = 0.75))#个数
p1 + p2# “gene+nt” 使用包含 TCR/Ig + CDR3 区域的核苷酸序列的 VDJC 基因
quantContig_output <- quantContig(combined, cloneCall="gene+nt", scale = T, 
                                  exportTable = T)#输出
quantContig_output
#2.克隆型的相对分布（丰度）
abundanceContig(combined, cloneCall = "gene", scale = F)+scale_fill_lancet()
#3.CDR3序列的长度分布，"nt"CDR3核苷酸序列"aa"氨基酸序列
lengthContig(combined, cloneCall="aa", chain = "single")+scale_fill_npg()
lengthContig(combined, cloneCall="nt")+scale_fill_npg()
#4.样本之间的克隆类型和动态变化，没有共享则没有连线
compareClonotypes(combined,numbers = 8,samples = c("MM1", "BB1"),
                  cloneCall="aa", graph = "alluvial") +  scale_fill_npg()
compareClonotypes(combined,numbers = 8,samples = c("MM1", "BB1"),
                  cloneCall="aa", graph = "alluvial") 
#5.克隆空间稳态
clonalHomeostasis(combined, cloneCall = "gene",
                  cloneTypes = c(Rare = 1e-04,
                                 Small = 0.001,
                                 Medium = 0.01,
                                 Large = 0.1,
                                 Hyperexpanded = 1))
#6.克隆比例，按克隆型的出现频率将其进行排名，1:10表示每个样本中的前10个克隆型
clonalProportion(combined, cloneCall = "nt")
#7. Overlap Analysis 重叠分析
clonalOverlap(combined, cloneCall = "gene+nt",   
              method = "overlap")

clonesizeDistribution(combined, 
                      cloneCall = "gene+nt", 
                      method="ward.D2")
#8.多样性分析
clonalDiversity(combined, 
                cloneCall = "gene+nt",  
                n.boots = 100)

clonalDiversity(combined, 
                cloneCall = "gene+nt",  
                group = "ID",
                n.boots = 100)
#和scRNA合并
scRNA <- combineExpression(combined, scRNA, 
                            cloneCall="gene", 
                            # group.by = "Sample", 
                            proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
saveRDS(seurat,'scRNA.TB6.rds')
DimPlot(scRNA, group.by = "cloneType",label = F)+
  DimPlot(scRNA, group.by = "seurat_clusters",label = T)
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTgene")] <- c("B.CTgene")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTnt")] <- c("B.CTnt")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTaa")] <- c("B.CTaa")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTstrict")] <- c("B.CTstrict")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("Frequency")] <- c("B.Frequency")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("cloneType")] <- c("B.cloneType")

colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTgene")] <- c("T.CTgene")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTnt")] <- c("T.CTnt")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTaa")] <- c("T.CTaa")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("CTstrict")] <- c("T.CTstrict")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("Frequency")] <- c("T.Frequency")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("cloneType")] <- c("T.cloneType")
colnames(scRNA@meta.data)[colnames(scRNA@meta.data) == c("barcode")] <- c("T.barcode")
