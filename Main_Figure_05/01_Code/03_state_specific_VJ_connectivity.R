Tcell$Celltype2 <- factor(Tcell$Celltype2, levels = c(
  "S1_CD4_Naive", "S1_CD4_Tscm", "S1_CD4_Tcm", "S1_CD4_Tem",  "S1_CD4_Treg",
  "S2_CD4_Naive", "S2_CD4_Tscm", "S2_CD4_Tcm", "S2_CD4_Tem", "S2_CD4_Teff", "S2_CD4_Treg", 
  "S3_CD4_Naive", "S3_CD4_Tscm", "S3_CD4_Tcm", "S3_CD4_Tem", "S3_CD4_Teff", "S3_CD4_Treg", 
  "S1_CD8_Naive", "S1_CD8_Tcm", "S1_CD8_Tem", "S1_CD8_Teff", "S1_CD8_Tex", 
  "S2_CD8_Naive", "S2_CD8_Tcm", "S2_CD8_Tem", "S2_CD8_Teff", "S2_CD8_Tex", 
  "S3_CD8_Naive", "S3_CD8_Tcm", "S3_CD8_Tem", "S3_CD8_Teff", "S3_CD8_Tex",
  "Trm-like", "NKT", "MAIT", "γδT"))

DimPlot(Tcell,group.by = 'Celltype2',cols = c(
  'white','white', 'white','white','white',
  'white','white', 'white','white','#A52A2A','white',
  'white','white', 'white','white','#4F7942','white',
  'white','white', 'white','white', 'white',
  'white','white', 'white', "#A0478D", 'white',
  'white','white', 'white', "#FFB300", 'white',
  'white','white', 'white','white'
))

compareClonotypes(CD8Tex,numbers=10,split.by = 'condition3',
                  samples=c("CD8_Teff_S2","Tex"),cloneCall="gene",graph="alluvial")


library(Seurat)
library(tidyverse)
library(ggalluvial)
library(RColorBrewer)

# 1. 提取并计算原始占比 (分母为该组总细胞数)
full_counts <- CD8Tex@meta.data %>%
  group_by(condition3) %>%
  summarise(total = n())

plot_df <- CD8Tex@meta.data %>%
  filter(condition3 %in% c("S2_CD8_Teff", "CD8_Tex")) %>%
  mutate(TRA_part = str_extract(CTgene, "TRAV[^.]+.[^.]+")) %>%
  filter(!is.na(TRA_part)) %>%
  group_by(condition3, TRA_part) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(full_counts, by = "condition3") %>%
  mutate(proportion = count / total)

# 2. 筛选 Top 10
top_10_vj <- plot_df %>%
  group_by(TRA_part) %>%
  summarise(sum_p = sum(proportion)) %>%
  slice_max(sum_p, n = 10) %>%
  pull(TRA_part)

final_df <- plot_df %>% filter(TRA_part %in% top_10_vj)

# 3. 排序与配色：确保从浅到深，且 Other 不干扰
vj_levels <- sort(top_10_vj)
# 将因子水平反转，使深色在下，浅色在上（匹配 ggplot 堆叠顺序）
final_df$TRA_part <- factor(final_df$TRA_part, levels = rev(vj_levels))
my_palette <- c("#6E5AB4", "#7D66BD", "#8C72C6", "#9B7ECF", "#AA8AD8", 
                "#B996E1", "#C8A2EA", "#D7AEF3", "#E6BAFC", "#F5C6FF")


# 4. 绘图：利用 Legend 模拟文献右侧注脚
ggplot(final_df, aes(x = condition3, y = proportion, 
                     stratum = TRA_part, alluvium = TRA_part)) +
  geom_stratum(aes(fill = TRA_part), color = "black", size = 0.3, width = 0.35) +
  geom_alluvium(aes(fill = TRA_part), alpha = 0.4, width = 0.35) +
  scale_fill_manual(values = rev(my_palette), name = "V-J Pairs") + # 这里设置图例标题
  theme_classic() +
  labs(x = "", y = "Proportion") +
  theme(
    # --- 关键：定制图例使其看起来像文献注脚 ---
    legend.position = "right",           # 放在右侧
    legend.title = element_blank(),      # 去掉图例标题，更干净
    legend.text = element_text(size = 10, face = "plain"), 
    legend.key.size = unit(0.5, "cm"),   # 控制色块大小
    legend.spacing.y = unit(0.2, "cm"),  # 增加标签之间的垂直间距，防止重叠
    
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.y = element_text(size = 10, color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  # 引导图例垂直排列
  guides(fill = guide_legend(reverse = FALSE, ncol = 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5))



library(Seurat)
library(tidyverse)
library(ggalluvial)

# 1. 提取并计算原始占比
full_counts <- CD8Tex@meta.data %>%
  group_by(condition3) %>%
  summarise(total = n())

plot_df <- CD8Tex@meta.data %>%
  filter(condition3 %in% c("S3_CD8_Teff", "CD8_Tex")) %>%
  # --- 关键修改区：提取并清洗字符串 ---
  mutate(TRB_part = str_extract(CTgene, "TRBV[0-9A-Z.-]+.TRBJ[0-9A-Z.-]+")) %>%
  filter(!is.na(TRB_part)) %>%
  mutate(TRB_part = TRB_part %>% 
           str_remove_all("\\.NA") %>%      # 移除 .NA
           str_remove_all("\\.TRBC[12]") %>% # 移除 .TRBC1 或 .TRBC2
           str_replace_all("\\.", ",")) %>%  # 将中间的点替换为逗号，匹配文献风格
  # ----------------------------------
group_by(condition3, TRB_part) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(full_counts, by = "condition3") %>%
  mutate(proportion = count / total)

# 2. 筛选 Top 10
top_10_vj <- plot_df %>%
  group_by(TRB_part) %>%
  summarise(sum_p = sum(proportion)) %>%
  slice_max(sum_p, n = 10) %>%
  pull(TRB_part)

final_df <- plot_df %>% filter(TRB_part %in% top_10_vj)

# 3. 排序与配色
vj_levels <- sort(top_10_vj)
final_df$TRB_part <- factor(final_df$TRB_part, levels = rev(vj_levels))

# 你要求的蓝绿渐变色调
my_palette <- c("#168D6D", "#1A9A78", "#1EA783", "#22B48E", "#26C199", 
                "#44CC9F", "#62D7A5", "#80E2AB", "#9EEDB1", "#BCF8B7")

# 4. 绘图
ggplot(final_df, aes(x = condition3, y = proportion, 
                     stratum = TRB_part, alluvium = TRB_part)) +
  geom_stratum(aes(fill = TRB_part), color = "black", size = 0.3, width = 0.35) +
  geom_alluvium(aes(fill = TRB_part), alpha = 0.4, width = 0.35) +
  scale_fill_manual(values = rev(my_palette)) + 
  theme_classic() +
  labs(x = "", y = "Proportion") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    axis.text.x = element_text(size = 10, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  ) +
  guides(fill = guide_legend(reverse = FALSE, ncol = 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0.5))