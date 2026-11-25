CD4@meta.data$TOX <- CD4[["RNA"]]@data["TOX", ]
# 计算每个 State2 中 TOX>0 的细胞占比
meta_df <- CD4@meta.data[, c("State2", "TOX")]
TOX_ratio <- tapply(meta_df$TOX > 0, meta_df$State2, mean)
plot_data <- data.frame(
  Group = names(TOX_ratio),
  Ratio = as.numeric(TOX_ratio)
)
plot_data$Group <- factor(plot_data$Group, levels = c("S1_CD4_Naive",'S1_CD4_Tscm',"S1_CD4_Tcm", "S1_CD4_Tem","S1_CD4_Treg", 
                                                      "S2_CD4_Naive",'S2_CD4_Tscm',"S2_CD4_Tcm", "S2_CD4_Tem", "S2_CD4_Teff", "S2_CD4_Treg", 
                                                      "S3_CD4_Naive",'S3_CD4_Tscm',"S3_CD4_Tcm", "S3_CD4_Tem", "S3_CD4_Teff", "S3_CD4_Treg"), ordered = TRUE)
custom_colors <- c('#4682B4', '#4682B4', '#4682B4', '#4682B4', '#4682B4', 
                   '#4682B4', '#4682B4', '#4682B4', '#4682B4','#A52A2A','#967969',
                   '#4682B4', '#4682B4', '#4682B4', '#4682B4','#4F7942','#097969')

# 绘制柱形图
ggplot(plot_data, aes(x = Group, y = Ratio)) +
  geom_col(aes(fill = Group), width = 0.7) +  # 使用aes映射fill到Group
  scale_fill_manual(values = custom_colors) +  # 设置自定义颜色
  geom_text(aes(label = percent(Ratio, accuracy = 0.1)), 
            vjust = -0.5, size = 2.5) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.05))
  ) +labs(
    title = "",x = "",y = "Ratio of TOX+ cell") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 8, face = "bold"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.3),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    # 添加实线XY轴
    axis.line = element_line(color = "black", linewidth = 0.3),
    # 去掉底框
    panel.border = element_blank(),
    # 确保轴线显示完整
    axis.ticks = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length = unit(0.2, "cm"),
    # 去除刻度线向图内延伸的灰色虚线
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    # 可选：隐藏图例
    legend.position = "none")
