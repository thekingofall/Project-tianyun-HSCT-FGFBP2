# 整体------------------------
DimPlot(scRNA,raster = F,label = T,group.by = 'Celltype',
        cols = c("#49B2BE",'#699ECA','#0098E2',"#FB5039",'#C07DD9','#D97DAA','#FF8C00',"#DE3163","#8582BD"))#,
DotPlot(scRNA,features = c("CD34","SPINK2","PROM1",
                           "CD3E","CD3D","CD3G",
                           "GNLY", "NKG7", "CX3CR1",
                           'CD79A','CD19','MS4A1',
                           "CD14","VCAN","LYZ",
                           "FCGR3B","CMTM2","S100A8",
                           "CD1C","FCER1A","ENHO",
                           'IL3RA','CLEC4C','JCHAIN',
                           "PPBP", "PF4", "GP9"),group.by = 'Celltype')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#密度图
p <- DimPlot(Tcell,raster = F,label = F,group.by = 'Celltype',split.by = 'condition',
             cols = c("#49B2BE",'#699ECA','#0098E2',"#FB5039",'#C07DD9','#D97DAA','#FF8C00',"#DE3163","#8582BD"))
p$data%>% 
  ggplot(aes(UMAP_1, UMAP_2)) +
  stat_density_2d(aes(fill=(..level..-min(..level..))/(max(..level..)-min(..level..))),
                  geom='polygon',contour=T) +
  scale_fill_distiller(palette="Blues",direction=1,name="Density") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),  # 设置分面图的外框为黑色
        strip.background = element_blank(),
        strip.text = element_text(size=12,hjust=0.5),
        legend.position="none") +
  facet_wrap(~condition, ncol = 5)
#细胞比例
mycolor<-c("#49B2BE",'#699ECA','#0098E2',"#FB5039",'#C07DD9','#D97DAA','#FF8C00',"#DE3163","#8582BD") 
t <- prop.table(table(scRNA$Celltype, scRNA$condition), margin = 2)
p <- as.data.frame(t)
#p <- p[!p$Var1 %in% c("Erythrocytes", "Platelets"), ]
p$x <- rep(c(1:5), each=9)
#p$Var1 <- factor(p$Var1, levels = 1:22)
ggplot(p, aes(x, Freq, fill = Var1)) +
  geom_area(size = 0.5, colour = "white") +
  scale_fill_manual(name = "Cluster", values = mycolor) +
  scale_x_continuous(breaks = 1:5, labels = c('HD','Pre','1M','3M','6M')) +
  theme_bw() +
  theme(panel.border = element_rect(size = 1),
        axis.text.x = element_text(size = 10),# 调整x轴文本大小为12,angle = 90,hjust = 0.5,vjust = 0.5
        axis.title.x = element_blank(),  # 去掉x轴标题
        axis.title.y = element_text(size = 12),  # 调整y轴标题字体大小为14
        axis.text.y = element_text(size = 10)) +  # 调整y轴文本字体大小为12
  labs(y = "Proportion (%)") +
  guides(fill = guide_legend(ncol = 1))

mycolor<-c('#7030a0','#0070c0','#ff5050','#f78832','#70ad47') 
t <- prop.table(table(scRNA$condition, scRNA$condition1), margin = 2)

mycolor <- c('#7030a0','#0070c0','#ff5050','#f78832','#70ad47')
t <- as.matrix(prop.table(table(scRNA$condition, scRNA$condition1), margin = 2)*100) 
p <- as.data.frame(t)
p$x <- rep(1:1, each=5)
colnames(p) <-c("cluster","Disease","proportion","x")
ggplot(p,aes(Disease,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =mycolor)+
  ylab("Proportion(%)")
#各患者占比
PS8 <- scRNA[,scRNA$orig.ident %in% c("GZ8", "YZ8", "YZ8A", "YZ8C", "YZ8D")]
PS7 <- scRNA[,scRNA$orig.ident %in% c("GZ7", "YZ7", "YZ7A", "YZ7C", "YZ7D")]
PS5 <- scRNA[,scRNA$orig.ident %in% c("GZ5", "YZ5", "YZ5A", "YZ5C", "YZ5D")]
PS6 <- scRNA[,scRNA$orig.ident %in% c("GZ6", "YZ6", "YZ6A", "YZ6C", "YZ6D")]
t <- as.matrix(prop.table(table(PS8$Celltype,PS8$condition), margin = 2)*100)#去除γδT
p <- as.data.frame(t)
colnames(p) <-c("cluster","Disease","proportion")
group.colors <- c("#49B2BE",'#699ECA','#0098E2',"#FB5039",'#C07DD9','#D97DAA','#FF8C00',"#DE3163","#8582BD")
ggplot(p,aes(Disease,proportion,fill=cluster))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5,vjust = 0.5))+
  guides(fill=guide_legend(title=NULL))+
  scale_fill_manual(values =group.colors)+ggtitle("Patient8")+
  labs(x = NULL, y = NULL)+
  theme(plot.title = element_text(hjust = 0.5))+coord_flip()
#血常规数据
# Load necessary libraries
library(ggplot2)
library(reshape2)
# Create the data frame
data <- data.frame(
  Patient = c("Patient5", "Patient6", "Patient7", "Patient8"),
  Pre_WBC = c(5.24, 5.85, 2.48, 5.5),
  Pre_LYMPH = c(0.72, 0.81, 0.69, 0.91),
  Post_1M_WBC = c(2.76, 4.39, 9.88, 7.92),
  Post_1M_LYMPH = c(0.6, 0.67, 0.73, 0.78),
  Post_3M_WBC = c(6.46, 3.27, 3.15, 6.0),
  Post_3M_LYMPH = c(1.56, 0.48, 0.5, 0.29),
  Post_6M_WBC = c(4.64, 2.33, 2.95, 3.78),
  Post_6M_LYMPH = c(1.08, 0.95, 1.02, 0.25)
)
# Melt the data for ggplot
data_melted <- melt(data, id.vars = "Patient")
# Create a new column for Time and Type
data_melted$Time <- rep(c("Pre", "Post-1M", "Post-3M", "Post-6M"), each = 4)
data_melted$Type <- rep(c("WBC", "LYMPH"), times = 8)
# Extract the numeric values for WBC and LYMPH
data_melted$Value <- NA
for (i in 1:nrow(data_melted)) {
  if (grepl("WBC", data_melted$variable[i])) {
    data_melted$Value[i] <- data_melted$value[i]
  } else {
    data_melted$Value[i] <- data_melted$value[i]
  }
}
# Specify colors for patients
patient_colors <- c("Patient5" = "#7030a0", "Patient6" = "#0070c0", 
                    "Patient7" = "#ff5050", "Patient8" = "#f78832")
ggplot(data_melted, aes(x = Time, y = Value, color = Patient)) +
  geom_point(aes(shape = Type), size = 3) +
  scale_shape_manual(values = c("WBC" = 16, "LYMPH" = 17)) +  # 16 is filled circle, 17 is triangle
  scale_color_manual(values = patient_colors) +
  labs(x = "Time", y = "Value") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Load necessary libraries
library(ggplot2)
library(reshape2)
# Create the data frame
data <- data.frame(
  Patient = c("Patient5", "Patient6", "Patient7", "Patient8"),
  Pre_WBC = c(5.24, 5.85, 2.48, 5.5),
  Pre_LYMPH = c(0.72, 0.81, 0.69, 0.91),
  Post_1M_WBC = c(2.76, 4.39, 9.88, 7.92),
  Post_1M_LYMPH = c(0.6, 0.67, 0.73, 0.78),
  Post_3M_WBC = c(6.46, 3.27, 3.15, 6.0),
  Post_3M_LYMPH = c(1.56, 0.48, 0.5, 0.29),
  Post_6M_WBC = c(4.64, 2.33, 2.95, 3.78),
  Post_6M_LYMPH = c(1.08, 0.95, 1.02, 0.25)
)
data_melted <- melt(data, id.vars = "Patient")
data_melted$Time <- factor(sub("_(WBC|LYMPH)", "", data_melted$variable), 
                           levels = c("Pre", "Post_1M", "Post_3M", "Post_6M"))
data_melted$Type <- sub(".*_", "", data_melted$variable)  # Extract WBC or LYMPH
colnames(data_melted)[3] <- "Value"
patient_colors <- c("Patient5" = "#621BAD", "Patient6" = "#1B3FAD",
                    "Patient7" = "#33AD1B", "Patient8" = "#D93119")
ggplot(data_melted, aes(x = Time, y = Value, color = Patient)) +
  geom_point(aes(shape = Type), size = 3) +
  scale_shape_manual(values = c("WBC" = 16, "LYMPH" = 17)) +  # 16 is filled circle, 17 is triangle
  scale_color_manual(values = patient_colors) +
  labs(x = "", y = "Peripheral Blood Cell Count") +
  theme_bw() +
  theme(axis.ticks.length = unit(0.1, 'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# T细胞------------------------
#克隆丰度
library(ggplot2)
library(dplyr)
p <- abundanceContig(Tcell,cloneCall = "CTgene",split.by = 'condition',exportTable = T)
hap <- Tcell[,Tcell$orig.ident %in% c("GZ5", "YZ5", "YZ5A", "YZ5C", "YZ5D",
                          "GZ6", "YZ6", "YZ6A", "YZ6C", "YZ6D", 
                          "GZ7", "YZ7", "YZ7A", "YZ7C", "YZ7D")]
P5 <- Tcell[,Tcell$orig.ident %in% c("GZ5", "YZ5", "YZ5A", "YZ5C", "YZ5D")]
P6 <- Tcell[,Tcell$orig.ident %in% c("GZ6", "YZ6", "YZ6A", "YZ6C", "YZ6D")]
P7 <- Tcell[,Tcell$orig.ident %in% c("GZ7", "YZ7", "YZ7A", "YZ7C", "YZ7D")]
P9 <- Tcell[,Tcell$orig.ident %in% c("GZ9", "YZ9", "YZ9A")]
mud <- Tcell[,Tcell$orig.ident %in% c("GZ8", "YZ8", "YZ8A", "YZ8C", "YZ8D")]

abundanceContig(P5,cloneCall = "CTaa",split.by = 'condition')+#,exportTable = T
  scale_color_manual(values = c('#0070c0','#ff5050','#f78832','#70ad47','#7030a0'))+
  ggtitle("P5 Clonal Abundance") +theme(plot.title = element_text(hjust = 0.5))
#克隆共享
P8 <- Tcell[,Tcell$condition1 %in% c("Patient8")]
P7 <- Tcell[,Tcell$condition1 %in% c("Patient7")]
P5 <- Tcell[,Tcell$condition1 %in% c("Patient5")]
P6 <- Tcell[,Tcell$condition1 %in% c("Patient6")]
clonalOverlap(P8,split.by = "condition",
              cloneCall="CTaa", method="overlap"#overlap raw
)+ scale_fill_gradient(low = "#3889E9", high = "#E9E638", na.value = "white")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5))+ggtitle("Single")+
  theme(plot.title = element_text(hjust = 0.5))+ggtitle("Patient8 MUD-HSCT") +  # 修改这里以添加您的标题
  theme(plot.title = element_text(hjust = 0.5)) #Haplo
abundanceContig(P8,cloneCall = "CTaa",split.by = 'condition')+#,exportTable = T
  scale_color_manual(values = c('#0070c0','#ff5050','#f78832','#70ad47','#7030a0'))+
  ggtitle("Patient8") +theme(plot.title = element_text(hjust = 0.5))

library(reticulate)
use_condaenv("cyto")  # 指向您创建的环境
#T细胞克隆多样性
library(ggplot2)
library(tidyr)
library(dplyr)
# 使用实际的数据
p <- clonalDiversity(Tcell,split.by = 'condition',
                     cloneCall = "CTaa",
                     n.boots = 100,exportTable = T)
# 重新排序Group
p$Group <- factor(p$Group, levels = c( "Pre","HD","Post-1M","Post-3M","Post-6M"))
# 将数据转换为长格式
p_long <- p %>%
  pivot_longer(cols = c(Shannon, Inv.Simpson, Chao, ACE, Inv.Pielou),
               names_to = "Index", 
               values_to = "Value")
# 确保Index顺序与原图一致
p_long$Index <- factor(p_long$Index, 
                       levels = c("Shannon", "Inv.Simpson", "Chao", "ACE", "Inv.Pielou"))
# 颜色映射
colors <- c("Pre" = '#7030a0',"HD" = '#0070c0', "Post-1M" = '#ff5050', 
            "Post-3M" = '#f78832', "Post-6M" = '#70ad47')
# 绘制连线图
ggplot(p_long, aes(x = Group, y = Value, color = Group)) +
  geom_point(size = 4) +
  geom_line(aes(group = 1), size = 1) +
  facet_wrap(~ Index, scales = "free_y", ncol = 5) +
  scale_color_manual(values = colors) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 10, face = "bold"),
    #panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次级网格线
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    legend.text = element_text(size = 8)
  ) +
  labs(
    title = "T Cell Clonal Diversity",
    y = "Index Score",
    x = ""
  )
tcr_list <- unique(unlist(data1))
Tcell$Tcb3 <- ifelse(Tcell$Tcb3 %in% tcr_list, Tcell$Tcb3, "other")
Tcell$Tcb4 <- ifelse(Tcell$Tcb3 %in% tcr_list, "WT1", "other")
#方块热图
library(ComplexHeatmap)
DefaultAssay(CD8T) <- "RNA"
features=c("PRKY", "UTY", "USP9Y", "IKZF2","ETFRF1","FASLG", "SPON2", "CERK", "S1PR5","SATB1","NEPRO",
           "XIST", "HLA-DQA1", "HLA-DRB5", "HLA-DRA", "HLA-DRB1", "HLA-DPB1", "HLA-DQB1", "GZMK", "AL138963.4","JUN","DUSP2")
gene_cell_exp <- AverageExpression(CD8T,
                                   features = features,
                                   group.by = 'condition',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'Group'
top_anno = HeatmapAnnotation(df = df,#细胞名/clusTcellr
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             gap = unit(5, "points"),
                             col = list(class = c("Pre"='#7030a0',
                                                  "Post-1M"='#ff5050',"Post-3M"='#f78832',"Post-6M"='#70ad47')))#颜色设置

marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
Heatmap(marker_exp,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(
          title=' '),
        col = colorRampPalette(c("navy","#E6E3A5","#E5E512"))(50), #'#008080',"white",'#CD5C5C'
        border = 'white',
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 3),
        top_annotation = top_anno,
        row_title_side = "right")
#T细胞功能评分
CD8CTL <- Tcell[,Tcell$Celltype %in% c("CD8 Tem","CD8 Teff","CD8 Tex")]
CD8TE <- Tcell[,Tcell$Celltype %in% c("CD8 Teff")]
CD8Na <- Tcell[,Tcell$Celltype %in% c("CD8 Naive")]
CD8Tex <- Tcell[,Tcell$Celltype %in% c("CD8 Tex")]
Palette <- c('#7030a0','#0070c0','#ff5050','#f78832','#70ad47')
my_comparisons <- list(c("Pre", "HD"),c('HD','Post-1M'),c('HD','Post-3M'),c('HD','Post-6M'))
plot_data <- data.frame(group = CD8Tex$condition, Aging = CD8Tex$texterm1)
ggplot(plot_data, aes(x = group, y = Aging, color = group, fill = group)) +
  geom_violin(position = position_dodge(width = 1), scale = 'width', color = "black") +
  geom_boxplot(position = position_dodge(width = 1), outlier.size = 0.7, width = 0.1, show.legend = FALSE, color = "black") +
  theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(size = 12),  # 设置x轴字体大小为12
        axis.text.y = element_text(size = 12),# 设置y轴字体大小为12
        plot.title = element_text(hjust = 0.5)) +  # 设置标题居中
  labs(x = '', y = '', title = 'CD8 Tex Exhaustion Score')+ #Cytotoxicity Signaling Activation
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",tip.length=0)+
  scale_fill_manual(values = Palette)

#柱状图
data <- data.frame(
  Patient = c("Patient5", "Patient6", "Patient7", "Patient8"),
  Pre = c(5366, 3284, 4349, 2684),
  HD = c(1518, 6577, 3673, 5858),
  Post_1M = c(1094, 3726, 520, 1810),
  Post_3M = c(7535, 3559, 732, 381),
  Post_6M = c(4362, 6537, 2347, 958)
)
#T细胞比例
data <- data.frame(
  Patient = c("Patient5", "Patient6", "Patient7", "Patient8"),
  Pre = c(76.41697522, 36.02062082, 51.68152109, 42.56263876),
  HD = c(49.44625407, 67.17393525, 45.25070839, 51.67607622),
  Post_1M = c(13.44806392, 53.8283733, 7.72889417, 17.41054252),
  Post_3M = c(70.06043701, 36.06241767, 9.42934433, 5.82925337),
  Post_6M = c(60.80289936, 64.57571866, 33.00056243, 20.87145969)
)

# 转换数据格式以适应 ggplot2
data_melted <- melt(data, id.vars = "Patient")

# 指定填充颜色
fill_colors <- c('#7030a0', '#0070c0', '#ff5050', '#f78832', '#70ad47')
fill_colors <- c('#B169DF', '#69AEDF', '#DF6969', '#DFB169', '#55B657')
# 绘制柱状图
ggplot(data_melted, aes(x = Patient, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = fill_colors) +
  labs(x = "", y = "Proportion of TCells(%)", title = "") +#Number of Tcells
  theme(panel.grid = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        legend.title = element_blank(), legend.key = element_blank(),
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))


# 加载必要的库
library(ggplot2)
library(reshape2)
# 创建数据框
data <- data.frame(
  Patient = c("Patient5", "Patient6", "Patient7", "Patient8"),
  Pre = c(5366, 3284, 4349, 2684),
  HD = c(1518, 6577, 3673, 5858),
  Post_1M = c(1094, 3726, 520, 1810),
  Post_3M = c(7535, 3559, 732, 381),
  Post_6M = c(4362, 6537, 2347, 958)
)
# 转换数据格式以适应 ggplot2
data_melted <- melt(data, id.vars = "Patient")
# 指定线条颜色
line_colors <- c("Patient5" = "#003366",   # 深蓝色
                 "Patient6" = "#00509E",   # 中蓝色
                 "Patient7" = "#66A3D2",   # 浅蓝色
                 "Patient8" = "#FF0000")   # 红色
# 绘制折线图
ggplot(data_melted, aes(x = variable, y = value, group = Patient, color = Patient)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = line_colors) +
  labs(x = "Time Point", y = "Value", title = "Patient Data Over Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))



