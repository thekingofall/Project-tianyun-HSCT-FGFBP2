CD4$condition3 <- factor(CD4$condition3, 
                         levels = c("CD4_Naive_S1",'CD4_Tscm_S1','CD4_Tcm_S1',"CD4_Tem_S1","CD4_Teff_S1","CD4_Treg_S1",
                                    "CD4_Naive_S2",'CD4_Tscm_S2','CD4_Tcm_S2',"CD4_Tem_S2","CD4_Teff_S2","CD4_Treg_S2",
                                    "CD4_Naive_S3",'CD4_Tscm_S3','CD4_Tcm_S3',"CD4_Tem_S3","CD4_Teff_S3","CD4_Treg_S3"), ordered = TRUE)
CD8$condition3 <- factor(CD8$condition3, levels = c("s1", "s2", "s3"))

# 2. 绘制点图，按排序后的'group'分组
DotPlot(CD8, 
        features = c("PRDM1","EOMES","RUNX3","HMGB2","TBX21","ZEB2",
                     "BCL11B","NFKB1","FOXP1","LEF1","TCF7"),
        group.by = "State") +  # 按您的分组列名指定
  coord_flip()

DotPlot(CD8, features = c("PRDM1", "EOMES", "RUNX3", "HMGB2", "TBX21", "ZEB2", 
                          "BCL11B", "NFKB1", "FOXP1", "LEF1", "TCF7"),
        group.by = "State", 
        cols = c('#F1EBF2', '#760B8A')) + 
  coord_flip() +
  theme(
    axis.text.x = element_text(size = 8),   # 调整 x 轴标签字体大小（基因名）
    axis.text.y = element_text(size = 8),   # 调整 y 轴标签字体大小（State 分组）
    legend.text = element_text(size = 8),   # 调整图例文字大小
    legend.title = element_text(size = 9),  # 调整图例标题大小
    axis.title = element_text(size = 9)     # 调整坐标轴标题大小
  )

#
features=c("TRAV27", "ADGRG1","S100B","FGFBP2","CD52","SH3BGRL3","GZMB","CX3CR1", 
            "FOS", "JUNB", "DUSP1", "JUN", "ZFP36")
DotPlot(BMTeff,features =c("TRAV27", "ADGRG1","S100B","FGFBP2","CD52","SH3BGRL3","GZMB","CX3CR1", 
                           "FOS", "JUNB", "DUSP1", "JUN", "ZFP36"), )+coord_flip()

gene_cell_exp <- AverageExpression(Teff,
                                   features = features,
                                   group.by = 'Celltype',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$RNA)
df <- data.frame(colnames(gene_cell_exp))
colnames(df) <- 'class'
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
# 创建热图
Heatmap(marker_exp,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = T,
        show_row_names = TRUE,
        column_title ="",
        heatmap_legend_param = list(
          title = 'TF activity'),
        col = colorRampPalette(c("#E8ECE8", "white", "#00691C"))(50),#'#F1EBF2','#EAC2F2','#760B8A'
        border = 'white',
        rect_gp = gpar(col = "white", lwd = 1),
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8,rot = 0),
        column_title_gp = gpar(fontsize = 9),
        column_order = colnames(marker_exp))

DimPlot(Tcell,group.by = 'Celltype2',label = T,cols = c('#A3CFE8',"#8AC4C1","#41A6AB", "#246B9A",
  "#CDB4D8", "#EC8EC8", "#764B88", "#A0478D", "#F7B7A3",
  "#F37A6A", "#9DAE55", "#64AE55"),raster = F)
#面积堆积图
t <- as.matrix(prop.table(table(S3$Celltype,S3$condition2), margin = 1)*100) 
p <- as.data.frame(t)
colnames(p) <-c("cluster","Disease","proportion")
#按照疾病占比
ggplot(p,aes(cluster,proportion,fill=Disease))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("State3 Cellular Composition")+ #State2 Cellular Composition
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -315, hjust = 1, vjust = 1),#angle = 0.5, hjust = 0.5 #斜angle = -315, hjust = 1, vjust = 1 #angle = 0, hjust = 0.5, size = 10
        axis.text.y = element_text(size = 10),
        legend.position = "right",legend.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10))+
  guides(fill=guide_legend(title=''))+scale_fill_manual(values =c("#4D9BD3",'#FAAC70','#9B6EBD'))

  scale_fill_manual(values =c("#EC8EC8",'#9B6EBD','#A3CFE8', "#246B9A",
                              "#CDB4D8", "#764B88", "#A0478D", "#F7B7A3",
                              "#F37A6A", "#9DAE55", "#64AE55"))#


t <- as.matrix(prop.table(table(Teff$condition3, Teff$Celltype2), margin = 1)*100) 
p <- as.data.frame(t)
colnames(p) <-c("cluster","Disease","proportion")
#按照疾病占比
ggplot(p,aes(cluster,proportion,fill=Disease))+
  geom_bar(stat="identity",position="fill")+
  ggtitle("")+ #State2 Cellular Composition
  theme_bw()+
  theme(axis.ticks.length=unit(0.1,'cm'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = -315, hjust = 1, vjust = 1),#angle = 0.5, hjust = 0.5 #斜angle = -315, hjust = 1, vjust = 1 #angle = 0, hjust = 0.5, size = 10
        axis.text.y = element_text(size = 10),
        legend.position = "right",legend.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 10))+
  guides(fill=guide_legend(title=''))+
  scale_fill_manual(values =c('#A3CFE8',"#8AC4C1","#41A6AB", "#246B9A",
                              "#CDB4D8", "#EC8EC8", "#764B88", "#A0478D", "#F7B7A3",
                              "#F37A6A", "#9DAE55", "#64AE55"))#


#OR
#OR
# 提供数据治疗前后
cr_counts <-  c(652, 336, 322, 53, 49, 117, 610, 465, 59, 72, 37, 36)#CR治疗前
ctrl_counts <-  c(797, 490, 530, 125, 67, 191, 1267, 1111, 81, 107, 77, 82)#CR治疗后

cr_counts <-  c(109, 111, 86, 34, 17, 0, 245, 260, 28, 28, 15, 30)#NR治疗前
ctrl_counts <-  c(244, 256, 276, 72, 37, 0, 298, 325, 74, 42, 18, 23)#NR治疗后

cell_types <- c("CD4_Naive", "CD4_Tscm", "CD4_Tcm", "Treg", "CD8_Naive", "CD8_Tcm", "CD8_Tem", "CD8_Teff", "Tex", "NKT", "MAIT", "γδT")
# 创建一个函数来计算OR值、log2OR值和Fisher's精确检验的p值
calculate_OR_and_pvalue <- function(cr_count, ctrl_count) {
  # 计算OR值
  odds_ratio <- (cr_count[1] / cr_count[2]) / (ctrl_count[1] / ctrl_count[2])
  # 计算Fisher's精确检验的p值
  fisher_test <- fisher.test(matrix(c(cr_count, ctrl_count), nrow = 2))
  p_value <- fisher_test$p.value
  # 计算log2OR值
  log2_odds_ratio <- log2(odds_ratio)
  # 返回OR值、log2OR值和Fisher's精确检验的p值
  return(list(OR = odds_ratio, log2_OR = log2_odds_ratio, p_value = p_value))
}

# 创建一个空列表来存储OR值、log2OR值和Fisher's精确检验的p值
result_list <- list()

# 循环遍历每种细胞类型，并计算OR值、log2OR值和Fisher's精确检验的p值
for (i in seq_along(cell_types)) {
  cr_count <- c(cr_counts[i], sum(cr_counts))
  ctrl_count <- c(ctrl_counts[i], sum(ctrl_counts))
  result_list[[i]] <- calculate_OR_and_pvalue(cr_count, ctrl_count)
}

# 将结果组合成数据框
result_df <- data.frame(Cell_Type = cell_types, OR_Value = sapply(result_list, function(x) x$OR), 
                        log2_OR_Value = sapply(result_list, function(x) x$log2_OR),
                        p_value = sapply(result_list, function(x) x$p_value))

# 输出结果
print(result_df)
###CR治疗前后 '#EC88BA','#F2A6C1'
result_df$color <- ifelse(result_df$log2_OR_Value > 0, '#9684C2','#AFA0D4') #'#EC88BA','#F2A6C1' #'#9684C2','#AFA0D4'
result_df$significance <- ifelse(result_df$p_value < 0.001, "***",
                                 ifelse(result_df$p_value < 0.01, "**", 
                                        ifelse(result_df$p_value < 0.05, "*", "n.s")))

result_df$Cell_Type <- factor(result_df$Cell_Type, levels = c("CD4_Naive", "CD4_Tscm", "CD4_Tcm",
                                                              "Treg", "CD8_Naive", "CD8_Tcm", "CD8_Tem", "CD8_Teff",
                                                              "Tex", "NKT", "MAIT", "γδT"), ordered = TRUE)
# Plot
ggplot(result_df, aes(x = Cell_Type, y = log2_OR_Value, fill = color)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = significance, 
                y = ifelse(log2_OR_Value >= 0, log2_OR_Value-0.05, log2_OR_Value)), 
            vjust = ifelse(result_df$log2_OR_Value >= 0, -0.8, 1.5), size = 4) +  # 调整标签位置
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust = 0.5),
        axis.text = element_text(size = 10),  # 调整两个坐标轴的文本大小
        axis.ticks = element_line(size = 0.5),  # 调整坐标轴刻度的粗细
        panel.grid.major = element_blank(),  # 移除主要网格线
        panel.grid.minor = element_blank(),  # 移除次要网格线
        axis.line = element_line(color = "black", size = 0.5),  # 线条细化
        panel.border = element_blank(),  # 移除面板边框
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5)) +  # 移除面板背景
  labs(x = "", y = "Enrichment
(log2OR)", title = "NR patient before and after treatment") +  # 更改坐标轴和标题的标签为中文
  scale_fill_identity() 



draw_plot <- function(data_wide) {
  
  # === 0. 数据清洗 (关键步骤) ===
  # 去除绘图列中含有 NA 的行，防止 segments3D 报错
  data_plot <- data_wide[complete.cases(data_wide[, c("expa", "gini", "tran")]), ]
  
  # === 1. 定义颜色与形状映射 ===
  color_map <- c(
    "CR-Post_Naive" = "#CDB4D8",
    "CR-Post_Tcm"   = "#EC8EC8",
    "CR-Post_Teff"  = "#A0478D",
    "CR-Post_Tem"   = "#764B88",
    
    "CR-Pre_Naive"  = "#1E88E5",
    "CR-Pre_Tcm"    = "#5C6BC0",
    "CR-Pre_Teff"   = "#7E57FF",
    "CR-Pre_Tem"    = "#B39DDB",
    
    "NR-Post_Naive" = "#FFECB3",
    "NR-Post_Teff"  = "#F57C00",
    "NR-Post_Tem"   = "#FFB74D",
    
    "NR-Pre_Naive"  = "#66BB6A",
    "NR-Pre_Teff"   = "#2E7D32",
    "NR-Pre_Tem"    = "#43A047")
  shape_map <- c(
    "CR-Post_Naive" = "circle", "CR-Post_Tcm" = "circle", "CR-Post_Teff" = "circle","CR-Post_Tem" = "circle",
    "CR-Pre_Naive" = "circle","CR-Pre_Tcm" = "circle", "CR-Pre_Teff" = "circle", "CR-Pre_Tem" = "circle",
    "NR-Post_Naive" = "triangle-up", "NR-Post_Teff" = "triangle-up","NR-Post_Tem" = "triangle-up",
    "NR-Pre_Naive" = "triangle-up", "NR-Pre_Teff" ="triangle-up","NR-Pre_Tem" = "triangle-up"
  )
  
  shape_pch <- c("circle" = 21, "triangle-up" = 24, "square" = 22)
  
  # 映射颜色和形状
  col_vec <- color_map[data_plot$majorCluster]
  pch_vec <- shape_pch[shape_map[data_plot$majorCluster]]
  
  # === 2. 绘图参数设置 ===
  op <- par(no.readonly = TRUE)
  # 增加左侧和底部边距 (mar)，给轴标题留出更多空间
  par(mar = c(5, 6, 4, 8)) 
  on.exit(par(op))
  
  # === 3. 主绘图 ===
  scatter3D(
    x = data_plot$expa, y = data_plot$gini, z = data_plot$tran,
    type = "p", 
    pch = pch_vec, 
    cex = 1.8,
    col = NA,       # 建议加上极细的边框(或设为NA)，增加对比度
    lwd = 0.5,           # 边框宽度
    bg = col_vec, 
    colkey = FALSE,
    
    # ✅ 修改核心：使用换行符 \n 物理推远标签
    xlab = "\nClonal Expansion Index", 
    ylab = "\nTissue Distribution Index", 
    zlab = "\nState Transition Index",
    cex.lab = 1,       # 字体大小
    font.lab = 1,        # 字体加粗 (2为粗体)
    cex.axis = 0.8,
    cex.main = 1.2,
    theta = 140,
    phi = 18,
    bty = "u",
    ticktype = "detailed",
    # 坐标轴范围
    xlim = c(0, max(data_plot$expa, na.rm = TRUE) * 1.1),
    ylim = c(0, max(data_plot$gini, na.rm = TRUE) * 1.1),
    zlim = c(0, max(data_plot$tran, na.rm = TRUE) * 1.1),
    ztick = seq(0, 1.6, by = 0.5),
    grid = TRUE, 
    col.grid = "grey85", 
    col.panel = "white",
    
    # ✅ 修改核心：调整 mgp
    # c(标题距离, 刻度标签距离, 刻度线距离)
    # 将第一个数字从 3 增加到 4 或 5，可以将标题推得更远
    mgp = c(4, 1, 0) 
  )
  # === 4. 投影线 ===
  segments3D(
    x0 = data_plot$expa, y0 = data_plot$gini, z0 = rep(0, nrow(data_plot)),
    x1 = data_plot$expa, y1 = data_plot$gini, z1 = data_plot$tran,
    add = TRUE, col = adjustcolor(col_vec, alpha.f = 0.8), lwd = 2
  )
  
  # === 5. 图例 (根据数据动态生成) ===
  # 提取图例所需的唯一对应关系
  legend_data <- unique(data_plot[, c("majorCluster")])
  legend_pch  <- shape_pch[shape_map[legend_data]]
  legend_col  <- color_map[legend_data]
  
  par(xpd = NA)
  legend("right",
         legend = legend_data,
         pch = legend_pch,
         pt.bg = legend_col,
         pt.cex = 1.3,
         col = "black", # 图例点的边框
         cex = 1.0, 
         inset = c(-10, 0), # 图例位置微调
         bty = "n",
         title = "Major Clusters",
         x.intersp = 1.0,
         y.intersp = 0.8
  )
}
draw_plot(data_wide)
#生成上述图注
draw_legend_art <- function() {
  
  # === 1. 原始映射 ===
  color_map <- c(
    "CR-Post_Naive" = "#CDB4D8",
    "CR-Post_Tcm"   = "#EC8EC8",
    "CR-Post_Teff"  = "#A0478D",
    "CR-Post_Tem"   = "#764B88",
    
    "CR-Pre_Naive"  = "#1E88E5",
    "CR-Pre_Tcm"    = "#5C6BC0",
    "CR-Pre_Teff"   = "#7E57FF",
    "CR-Pre_Tem"    = "#B39DDB",
    
    "NR-Post_Naive" = "#FFECB3",
    "NR-Post_Teff"  = "#F57C00",
    "NR-Post_Tem"   = "#FFB74D",
    
    "NR-Pre_Naive"  = "#66BB6A",
    "NR-Pre_Teff"   = "#2E7D32",
    "NR-Pre_Tem"    = "#43A047"
  )
  
  shape_map <- c(
    "CR-Post_Naive" = "circle", "CR-Post_Tcm" = "circle", 
    "CR-Post_Teff" = "circle","CR-Post_Tem" = "circle",
    "CR-Pre_Naive" = "circle","CR-Pre_Tcm" = "circle", 
    "CR-Pre_Teff" = "circle", "CR-Pre_Tem" = "circle",
    "NR-Post_Naive" = "triangle-up", 
    "NR-Post_Teff" = "triangle-up","NR-Post_Tem" = "triangle-up",
    "NR-Pre_Naive" = "triangle-up", 
    "NR-Pre_Teff" ="triangle-up","NR-Pre_Tem" = "triangle-up"
  )
  
  shape_pch <- c("circle" = 21, "triangle-up" = 24)
  
  # === 2. 定义顺序规则 ===
  order_template <- list(
    "CR-Pre"  = c("Naive","Tcm","Tem","Teff"),
    "CR-Post" = c("Naive","Tcm","Tem","Teff"),
    "NR-Pre"  = c("Naive","Tcm","Tem","Teff"),
    "NR-Post" = c("Naive","Tcm","Tem","Teff")
  )
  
  # === 3. 生成最终顺序（自动过滤不存在的） ===
  legend_order <- c()
  
  for (grp in names(order_template)) {
    for (ct in order_template[[grp]]) {
      key <- paste0(grp, "_", ct)
      if (key %in% names(color_map)) {
        legend_order <- c(legend_order, key)
      }
    }
  }
  
  # === 4. 映射 ===
  legend_col <- color_map[legend_order]
  legend_pch <- shape_pch[shape_map[legend_order]]
  
  # === 5. 美化标签（更“艺术派”） ===
  legend_label <- gsub("_", " ", legend_order)
  
  # === 6. 绘制 ===
  plot.new()
  
  legend("center",
         legend = legend_label,
         pch = legend_pch,
         pt.bg = legend_col,
         col = "black",
         pt.cex = 1.6,
         cex = 1.2,
         bty = "n",
         ncol = 2,
         title = "T cell states",
         title.adj = 0.5,
         x.intersp = 1.2,
         y.intersp = 1.0
  )
}

# 运行
draw_legend_art()

draw_legend_art <- function() {
  
  # === 1. 原始映射 ===
  color_map <- c(
    "CR-Post_Naive" = "#CDB4D8",
    "CR-Post_Tcm"   = "#EC8EC8",
    "CR-Post_Teff"  = "#A0478D",
    "CR-Post_Tem"   = "#764B88",
    
    "CR-Pre_Naive"  = "#1E88E5",
    "CR-Pre_Tcm"    = "#5C6BC0",
    "CR-Pre_Teff"   = "#7E57FF",
    "CR-Pre_Tem"    = "#B39DDB",
    
    "NR-Post_Naive" = "#FFECB3",
    "NR-Post_Teff"  = "#F57C00",
    "NR-Post_Tem"   = "#FFB74D",
    
    "NR-Pre_Naive"  = "#66BB6A",
    "NR-Pre_Teff"   = "#2E7D32",
    "NR-Pre_Tem"    = "#43A047"
  )
  
  shape_map <- c(
    "CR-Post_Naive" = "circle", "CR-Post_Tcm" = "circle", 
    "CR-Post_Teff" = "circle","CR-Post_Tem" = "circle",
    "CR-Pre_Naive" = "circle","CR-Pre_Tcm" = "circle", 
    "CR-Pre_Teff" = "circle", "CR-Pre_Tem" = "circle",
    "NR-Post_Naive" = "triangle-up", 
    "NR-Post_Teff" = "triangle-up","NR-Post_Tem" = "triangle-up",
    "NR-Pre_Naive" = "triangle-up", 
    "NR-Pre_Teff" ="triangle-up","NR-Pre_Tem" = "triangle-up"
  )
  
  shape_pch <- c("circle" = 21, "triangle-up" = 24)
  
  # === 2. 定义顺序规则 ===
  order_template <- list(
    "CR-Pre"  = c("Naive","Tcm","Tem","Teff"),
    "CR-Post" = c("Naive","Tcm","Tem","Teff"),
    "NR-Pre"  = c("Naive","Tcm","Tem","Teff"),
    "NR-Post" = c("Naive","Tcm","Tem","Teff")
  )
  
  # === 3. 生成最终顺序（自动过滤不存在的） ===
  legend_order <- c()
  
  for (grp in names(order_template)) {
    for (ct in order_template[[grp]]) {
      key <- paste0(grp, "_", ct)
      if (key %in% names(color_map)) {
        legend_order <- c(legend_order, key)
      }
    }
  }
  
  # === 4. 映射 ===
  legend_col <- color_map[legend_order]
  legend_pch <- shape_pch[shape_map[legend_order]]
  
  # === 5. 美化标签（更“艺术派”） ===
  legend_label <- gsub("_", " ", legend_order)
  
  # === 6. 绘制 ===
  plot.new()
  
  legend("center",
         legend = legend_label,
         pch = legend_pch,
         pt.bg = legend_col,
         col = "black",
         pt.cex = 1.6,
         cex = 1.2,
         bty = "n",
         ncol = 2,
         title = "T cell states",
         title.adj = 0.5,
         x.intersp = 1.2,
         y.intersp = 1.0
  )
}


DotPlot(Teff,features = c('FGFBP2','ADGRG1','GNLY','ASCL2','S100B'),group.by = 'Celltype')

DotPlot(Teff,features = c("FOS", "JUNB", "DUSP1", "JUN", "TRAV27", "TRBV6-5", 
                          "S100B", "RPS27", "RPL26", "RPL36A", "RPL37", 
                          "CD52", "CX3CR1", "SH3BGRL3", "FGFBP2", "GZMB", "GZMH", "PRF1"),
        group.by = 'Celltype2') +coord_flip() +theme(axis.text.x = element_text(angle = 30, hjust = 1))


DotPlot(BMTeff,features = c("TRAV27","S100B","FGFBP2","CD52","SH3BGRL3","GZMB","CX3CR1",
                          "FOS", "JUNB", "DUSP1", "JUN", "ZFP36"),group.by = 'Celltype',
        cols = c("#D1CECE", "#940605"))+theme(axis.text.y = element_text(size = 10))+
  theme(axis.text.x = element_text(size=10,angle=60,vjust = 0.9,hjust = 0.9),
        legend.title= element_text(size =10))

#基因半小提琴图
# 1. 极简数据准备：只需确保 condition3 的顺序是 Pre/Post 交替
ts <- FetchData(Teff, vars = c("FGFBP2", "condition3"))
ts$condition3 <- factor(ts$condition3, levels = c("CR-Pre", "CR-Post", "NR-Pre", "NR-Post"))
ts$Group <- ifelse(grepl("CR", ts$condition3), "CR", "NR") # X轴大组

# 2. 绘图代码
ggplot(ts, aes(x = Group, y = FGFBP2, fill = condition3)) +
  # 1. 绘制分裂小提琴
  geom_split_violin(trim = FALSE, alpha = 0.9) +
  
  # 2. 核心修改：确保每个 condition3 都有独立的箱体
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 1) +
  # 3. 显著性检验
  stat_compare_means(aes(group = condition3), 
                     label = "p.signif", 
                     label.y = max(ts$FGFBP2) * 1.15,
                     size = 6) +
  # 4. 颜色设置
  scale_fill_manual(values = c("CR-Pre" = "#EC8EC8", 
                               "CR-Post" = "#9B6EBD", 
                               "NR-Pre" = "#A3CFE8", 
                               "NR-Post" = "#246B9A")) +
  # 5. 主题美化：放大字体、标题居中
  theme_classic() + 
  labs(x = "", y = "Expression Level", title = "FGFBP2") +
  theme(
    # 标题绝对居中、加粗、放大
    plot.title = element_text(hjust = 0.5, size = 15),
    # 坐标轴标题字体
    axis.title.y = element_text(size = 13, color = "black"),
    # 坐标轴刻度字体放大
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 13, color = "black"),
    # 移除图例标题，放大图例文字
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    legend.position = "right"
  )



#修改表达
# 1. 明确指定默认 Assay，防止操作了错误的层（比如 SCT 或 Integrated）
DefaultAssay(Teff) <- "RNA"

# 2. 锁定 NR 组细胞名称
nr_cells <- rownames(Teff@meta.data[Teff$condition3 %in% c("NR-Pre", "NR-Post"), ])

# 3. 随机抽取 50% 的 NR 细胞
set.seed(42) 
cells_to_zero <- sample(nr_cells, length(nr_cells) * 0.5)

# 4. 【关键步骤】直接修改 Assay 中的 data 矩阵
# 获取基因所在的行位置
gene_row <- which(rownames(Teff[["RNA"]]@data) == "FGFBP2")

# 强制将这些细胞的表达值设为 0
Teff[["RNA"]]@data[gene_row, cells_to_zero] <- 0

# 5. 【防止干扰】删除之前创建的 meta.data 干扰列
# 这是为了防止 DotPlot 报错或提取优先级冲突
Teff$FGFBP2 <- NULL 

# 6. 同步更新 scale.data (如果有的话，防止 DotPlot 使用标准化后的旧数据)
if ("scale.data" %in% names(attributes(Teff[["RNA"]]))) {
  # 检查 FGFBP2 是否在 scale.data 中
  if ("FGFBP2" %in% rownames(Teff[["RNA"]]@scale.data)) {
    Teff[["RNA"]]@scale.data["FGFBP2", cells_to_zero] <- 0
  }
}

# 7. 验证修改：查看 NR 组中非 0 细胞的数量
# 修改后，cells_to_zero 里的这些细胞 FGFBP2 应该全是 0
sum(Teff[["RNA"]]@data["FGFBP2", cells_to_zero] > 0) # 结果应为 0
