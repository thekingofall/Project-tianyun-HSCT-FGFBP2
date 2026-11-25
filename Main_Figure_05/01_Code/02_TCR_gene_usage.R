# 确保是字符型
ctgene <- as.character(P5$CTgene)

# 先按 α / β 链拆分
split_ab <- strsplit(ctgene, "_")

# 初始化结果向量
TRAV <- TRAJ <- TRBV <- TRBJ <- rep(NA_character_, length(split_ab))

for (i in seq_along(split_ab)) {
  
  parts <- split_ab[[i]]
  
  # α链
  if (length(parts) >= 1 && !is.na(parts[1])) {
    alpha <- strsplit(parts[1], "\\.")[[1]]
    TRAV[i] <- ifelse(length(alpha) >= 1, alpha[1], NA)
    TRAJ[i] <- ifelse(length(alpha) >= 2, alpha[2], NA)
  }
  
  # β链
  if (length(parts) >= 2 && !is.na(parts[2])) {
    beta <- strsplit(parts[2], "\\.")[[1]]
    TRBV[i] <- ifelse(length(beta) >= 1, beta[1], NA)
    TRBJ[i] <- ifelse(length(beta) >= 3, beta[3], NA)
  }
}

# 写入 Seurat metadata
P5$TRAV <- TRAV
P5$TRAJ <- TRAJ
P5$TRBV <- TRBV
P5$TRBJ <- TRBJ

df <- table(P5$TRBV,P5$condition3)
write.csv(df,'./11月富集/P5TRBV.csv')

#绘制靶心玫瑰图
library(ggseqlogo)
library(ggplot2)
library(dplyr)

df <-  data.frame(
  TRAJ = c("TRAJ12", "TRAJ4", "TRAJ30", "TRAJ45", "TRAJ9", "TRAJ38", "TRAJ17", "TRAJ49", "TRAJ32", "TRAJ33"),
  freq = c(27, 19, 13, 11, 11, 8, 7, 7, 7, 7),
  prop = c(23.1, 16.2, 11.1, 9.4, 9.4, 6.8, 6.0, 6.0, 6.0, 6.0)
)

# 核心修改：按 freq 从小到大排序，这样玫瑰花瓣会由小到大排列
df <- df %>% arrange(freq)
# 锁定因子顺序，防止 ggplot 自动按字母重排
df$TRAJ <- factor(df$TRAJ, levels = df$TRAJ)

# 2. 绘图
ggplot(df, aes(x = TRAJ, y = freq, fill = TRAJ)) +
  # 设置宽度，让扇形之间有微小间隙或紧凑排列
  geom_col(width = 0.9, alpha = 0.8, color = "white") + 
  
  # coord_polar 核心设置：
  # start=0 控制起始角度
  coord_polar(start = 0) +
  
  # 1. 实现中心空心：通过 ylim 的下限设置为负数来实现
  # 负值越大，中间空心越大。这里设为 -500 左右较为适中
  ylim(-100, max(df$freq) + 200) + 
  
  # 2. 外部标签与引线：使用 geom_text_repel 或 geom_label_repel
  # 也可以直接用 geom_text 手动计算位置
  geom_text_repel(
    aes(y = freq, label = paste0(TRAJ, "\n", prop, "%")),
    nudge_y = 150,           # 标签向外推的距离
    segment.size = 0.3,      # 引线的粗细
    segment.color = "#E3E6E6", # 引线的颜色
    direction = "y",         # 强制垂直方向排布防止重叠
    size = 3#fontface = "bold"
  ) +
  
  # 3. 颜色与主题
  scale_fill_brewer(palette = "Set3") +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = margin(1, 1, 1, 1, "cm") # 留白防止标签切断
  )

