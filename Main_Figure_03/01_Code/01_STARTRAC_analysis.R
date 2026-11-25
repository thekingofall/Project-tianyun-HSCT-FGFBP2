#head(subT.meta)
#Clone_NUM orig.ident condition condition1 condition2  Celltype State majorCluster stype
#1         1       YZ8A   Post-1M   Patient8       Post  CD8_Teff    S3  S3_CD8_Teff   CD8
#2         1        GZ8        HD   Patient8         HD CD8_Naive    S1 S1_CD8_Naive   CD8
#3         1       YZ5A   Post-1M   Patient5       Post  CD8_Teff    S1  S1_CD8_Teff   CD8
#4         1        GZ6        HD   Patient6         HD   CD8_Tcm    S1   S1_CD8_Tcm   CD8
#5         1       YZ5D   Post-6M   Patient5       Post  CD8_Teff    S1  S1_CD8_Teff   CD8
#6         1       YZ5C   Post-3M   Patient5       Post   CD8_Tcm    S3   S3_CD8_Tcm   CD8
#clone.status patient sampleType               Cell_Name clone.id  loc
#1     NoClonal    YZ8A       Post YZ8A_CTGGTCTCACCATGTA-1  Clone_1 PBMC
#2     NoClonal     GZ8         HD  GZ8_GGGATGAGTAGTACCT-1  Clone_2 PBMC
#3     NoClonal    YZ5A       Post YZ5A_GAACCTAGTTAAGTAG-1  Clone_3 PBMC
#4     NoClonal     GZ6         HD  GZ6_TACACGAAGCAGATCG-1  Clone_4 PBMC
#5     NoClonal    YZ5D       Post YZ5D_AACTTTCGTAACGCGA-1  Clone_5 PBMC
#6     NoClonal    YZ5C       Post YZ5C_GAAGCAGAGTCTTGCA-1  Clone_6 PBMC
#合成一个新的
sample_names <- c(
  "AML1002.1.GEX","AML1002.4.GEX", 
  "AML1008.1.GEX","AML1008.3.GEX", 
  "AML1010.1.GEX","AML1010.3.GEX", 
  "AML1012.1.GEX","AML1012.3.GEX",
  "AML1016.1.GEX","AML1016.2.GEX",
  "AML1019.1.GEX","AML1019.3.GEX",
  "AML1022.1.GEX","AML1022.3.GEX",
  "AML1026.1.GEX","AML1026.3.GEX"
)

# 给每个data.frame加barcode前缀 + sample列
contig_list2 <- lapply(seq_along(contig_list), function(i) {
  df <- contig_list[[i]]
  
  df$barcode <- paste0(sample_names[i], "_", df$barcode)
  df$sample  <- sample_names[i]
  
  return(df)
})

# 合并（按列名自动对齐，缺失补NA）
vdj <- bind_rows(contig_list2)


vdj <- vdj %>%   
  dplyr::filter(high_confidence =="true" &                   
                                 chain %in% c("TRA","TRB") & 
                                 productive =="true")
vdj$Cell_name <- vdj$barcode
subT <- CD8
subT@meta.data <- subT@meta.data %>%   
  mutate(Cell_name = rownames(subT@meta.data)) %>%   
  inner_join(final_vdj, by = "Cell_name")
#结合TRA 和TRB的cdr3序列构建clone ，并统计每种clone的个数
subT@meta.data$Clone_AA = paste(subT@meta.data$cdr3.TRA, subT@meta.data$cdr3.TRB, sep="_")
subT@meta.data = subset(subT@meta.data, productive.TRA == "true" & productive.TRB == "true"  ) ; dim(subT@meta.data)
subT@meta.data = subT@meta.data %>% arrange(., Clone_AA)
### calculate clone number and clone ID
tmp = subT@meta.data %>%   
  group_by(Clone_AA) %>%  
  summarize(Clone_NUM = n()) %>%  
  mutate(Clone_ID = paste0("Clone_",rownames(.)))
head(tmp)
subT@meta.data = merge.data.frame(subT@meta.data, tmp) 

head(subT@meta.data,2)
subT@meta.data$clone.status <- ifelse(subT@meta.data$Clone_NUM >1 ,"Clonal","NoClonal")
subT@meta.data$Celltype <- gsub("CD4 ", "CD4_", subT@meta.data$Celltype)
subT@meta.data$majorCluster <- paste0(subT@meta.data$State,'_',subT@meta.data$Celltype)
subT@meta.data$stype <- 'CD4'
subT@meta.data$patient <- subT@meta.data$orig.ident
subT@meta.data$sampleType <- subT@meta.data$condition2
subT@meta.data$loc <- 'PBMC'
subT@meta.data$Cell_Name <- subT@meta.data$barcode
subT@meta.data$clone.id <- subT@meta.data$Clone_ID
subT.meta1 <- subT@meta.data %>%
  select(Cell_Name, clone.id,clone.status,patient,sampleType,stype,majorCluster,loc)
# 转换所有因子列为字符
subT.meta1 <- as.data.frame(lapply(subT.meta1, function(x) {
  if(is.factor(x)) {
    return(as.character(x))
  } else {
    return(x)
  }
}), stringsAsFactors = FALSE)
#计算
out2 <- Startrac.run(subT.meta1, proj="CRC",verbose=F)
p <- ggboxplot(as.data.table(out2@cluster.sig.data)[,][order(majorCluster),],
          x="majorCluster",y="value",palette = "npg",
          color = "index", add = "point", outlier.colour=NULL) +
  facet_wrap(~index,ncol=1,scales = "free_y") +
  theme(axis.text.x=element_text(angle = 60,hjust = 1))

#绘制3D散点图
library(plot3D)
library(RColorBrewer)

data_wide <- out2@cluster.data 
data_wide <- out2@cluster.data %>%
  group_by(majorCluster) %>%
  summarise(
    expa = median(expa, na.rm = TRUE),
    gini = median(gini, na.rm = TRUE),
    tran = median(tran, na.rm = TRUE),
    n_samples = n()  # 可选：添加样本数量
  ) %>%
  ungroup()

#本来配色
draw_plot <- function(data_wide) {
  # === 定义颜色与形状映射 ===
  color_map <- c( 
    "S1_CD4_Naive"='#B6D0E2', "S1_CD4_Tscm"='#6495ED', "S1_CD4_Tcm"='#0096FF',
    "S1_CD4_Tem"='#7393B3', "S1_CD4_Treg"='#5D3FD3',
    
    "S2_CD4_Naive"='#E5AA70', "S2_CD4_Tscm"='#C19A6B', "S2_CD4_Tcm"='#C4A484', 
    "S2_CD4_Tem"='#C04000',"S2_CD4_Teff"='#A52A2A', "S2_CD4_Treg"='#967969', 
    
    "S3_CD4_Naive"='#C1E1C1', "S3_CD4_Tscm"='#98FB98', "S3_CD4_Tcm"='#C9CC3F',
    "S3_CD4_Tem"='#4CBB17', "S3_CD4_Teff"='#4F7942', "S3_CD4_Treg"='#097969'
  )
  
  shape_map <- c(
    "S1_CD4_Naive"="circle", "S1_CD4_Tscm"="circle", "S1_CD4_Tcm"="circle", 
    "S1_CD4_Tem"="circle", "S1_CD4_Treg"="circle",
    
    "S2_CD4_Naive"="triangle-up", "S2_CD4_Tscm"="triangle-up", "S2_CD4_Tcm"="triangle-up", 
    "S2_CD4_Tem"="triangle-up", "S2_CD4_Treg"="triangle-up", "S2_CD4_Teff"="triangle-up",
    
    "S3_CD4_Naive"="square", "S3_CD4_Tscm"="square", "S3_CD4_Tcm"="square",
    "S3_CD4_Tem"="square", "S3_CD4_Teff"="square", "S3_CD4_Treg"="square"
  )
  shape_pch <- c("circle" = 21, "triangle-up" = 24, "square" = 22)
  
  clusters <- unique(data_wide$majorCluster)
  col_vec <- color_map[data_wide$majorCluster]
  pch_vec <- shape_pch[shape_map[data_wide$majorCluster]]
  
  # === 绘图参数 ===
  op <- par(no.readonly = TRUE)
  par(mar = c(5, 5, 2, 6))
  par(mgp = c(3, 1.5, 1.5))  # ✅ 增大轴标题间距
  par(tcl = 0)         # ✅ 缩短刻度线
  on.exit(par(op))
  
  # === 主绘图 ===
  scatter3D(
    x = data_wide$expa, y = data_wide$gini, z = data_wide$tran,
    type = "p", pch = pch_vec, cex = 1.8,col = NA,              
    bg = col_vec, colkey = FALSE,
    xlab = "Clonal Expansion Index", 
    ylab = "Tissue Distribution Index", 
    zlab = "State Transition Index",
    cex.lab = 1,     # ✅ 坐标轴标题字体更大
    font.lab = 1,      # ✅ 加粗
    cex.axis = 0.8,    # ✅ 坐标轴刻度字体更大
    cex.main = 1.2,
    theta = 140,
    phi = 18,
    bty = "u",
    ticktype = "detailed",
    xlim = c(0, max(data_wide$expa, na.rm = TRUE) * 1.1),
    ylim = c(0, max(data_wide$gini, na.rm = TRUE) * 1.1),
    zlim = c(0, 0.6),ztick = seq(0, 0.6, by = 0.2),
    lwd = 1.2, grid = TRUE, col.grid = "grey85", col.panel = "white",
    mgp = c(3, 3, 3) 
  )
  
  # === 投影线 ===
  segments3D(
    x0 = data_wide$expa, y0 = data_wide$gini, z0 = rep(0, nrow(data_wide)),
    x1 = data_wide$expa, y1 = data_wide$gini, z1 = data_wide$tran,
    add = TRUE, col = adjustcolor(col_vec, alpha.f = 0.8), lwd = 2
  )
  
  # === 图例 ===
  par(xpd = NA)
  legend("right",
         legend = clusters,
         pch = shape_pch[shape_map[clusters]],
         pt.bg = color_map[clusters],
         pt.cex = 1.3,
         col = NA,
         cex = 1.0, inset = c(0, 0), bty = "n",
         title = "Major Clusters",x.intersp = 1.8
  )
}
# 调用绘图
draw_plot(data_wide)

draw_plot <- function(data_wide) {
  
  # === 0. 数据清洗 (关键步骤) ===
  # 去除绘图列中含有 NA 的行，防止 segments3D 报错
  data_plot <- data_wide[complete.cases(data_wide[, c("expa", "gini", "tran")]), ]
  
  # === 1. 定义颜色与形状映射 ===
  color_map <- c( 
    "S1_CD8_Naive" = "#26A69A", "S1_CD8_Tcm" = "#00897B", "S1_CD8_Tem" = "#0097A7",
    "S1_CD8_Teff" = "#0288D1", "S1_CD8_Tex" = "#1565C0",
    "S2_CD8_Naive" = "#CDB4D8", "S2_CD8_Tcm" = "#EC8EC8", "S2_CD8_Tem" = "#764B88",
    "S2_CD8_Teff" = "#A0478D", "S2_CD8_Tex" = "#F7B7A3",
    "S3_CD8_Naive" = "#FFECB3", "S3_CD8_Tcm" = "#FFD54F", "S3_CD8_Tem" = "#FFCA28",
    "S3_CD8_Teff" = "#FFB300", "S3_CD8_Tex" = "#FF8F00"
  )
  
  shape_map <- c(
    "S1_CD8_Naive" = "circle", "S1_CD8_Tcm" = "circle", "S1_CD8_Tem" = "circle",
    "S1_CD8_Teff" = "circle", "S1_CD8_Tex" = "circle",
    "S2_CD8_Naive" = "triangle-up", "S2_CD8_Tcm" = "triangle-up", "S2_CD8_Tem" = "triangle-up",
    "S2_CD8_Teff" ="triangle-up","S2_CD8_Tex" = "triangle-up",
    "S3_CD8_Naive" = "square", "S3_CD8_Tcm" = "square", "S3_CD8_Tem" = "square",
    "S3_CD8_Teff" = "square", "S3_CD8_Tex" = "square"
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
    zlim = c(0, 1.6),
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

#--------------------------
p <- plot_ly(
  data = data_wide,x = ~expa, y = ~gini, z = ~tran,
  color = ~majorCluster,colors = color_map,
  symbol = ~majorCluster,symbols = shape_map,
  text = ~paste(
    "<b>", majorCluster, "</b>",
    "<br>扩增指数:", round(expa, 4),
    "<br>Gini指数:", round(gini, 4),
    "<br>转化指数:", round(tran, 4),
    "<br>样本数:", n_samples
  ),
  hoverinfo = "text",
  marker = list(size = 5, opacity = 1, line = list(width = 1, color = 'white'))
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(title = "Clonal Expansion Index", titlefont = list(size = 10)),
      yaxis = list(title = "Tissue Distribution Index", titlefont = list(size = 10)),
      zaxis = list(title = "State Transition Index", titlefont = list(size = 10)),
      camera = list(eye = list(x = 0.5, y = 0.5, z = 0.5))
    ),
    title = list(text = "STARTRAC Indices", x = 0.5, font = list(size = 14))
  )

p

