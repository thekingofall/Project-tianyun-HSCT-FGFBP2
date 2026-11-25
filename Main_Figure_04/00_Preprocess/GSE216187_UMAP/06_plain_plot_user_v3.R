library(Seurat)
library(ggplot2)

cell_colors <- c(
  "CD8 NV"="#DFC4DE", "CD8 hobit"="#b092b2",
  "CD8 mem. 1"="#d475a8", "CD8 mem. 2"="#a8c5db",
  "CD8 mem. 3"="#7eb5d6", "CD8 eff. 1"="#569ecf",
  "CD8 eff. 2"="#3d7ab8", "CD8 IFN"="#6b429e",
  "CD4 NV"="#e38888", "CD4 mem."="#e36442",
  "Th17"="#d9382b", "Treg"="#a81c1c", "CD4 IFN"="#7a1122",
  "γδT"="#8b52a3", "MAIT"="#5c2d82",
  "MLP"="#1f8a8a", "MEP"="#229eab", "NP"="#3bc4be",
  "MDP"="#42b883", "pDCs"="#7bc993", "cDCs"="#4db86e",
  "preB"="#2e9e52", "pro/pre T"="#1b7d3a", "MP/mono"="#115c27"
)

umap <- Embeddings(Bseurat, "umap")
df <- data.frame(
  UMAP_1=umap[, 1], UMAP_2=umap[, 2],
  Celltype=as.character(Bseurat@meta.data[rownames(umap), "Celltype"]),
  row.names=rownames(umap)
)
df$Celltype[df$Celltype %in% c("gdT", "γδ T")] <- "γδT"
df$Celltype <- factor(df$Celltype, levels=names(cell_colors))

# Label positions for this UMAP; the order follows cell_colors.
label_pos <- data.frame(
  Celltype=names(cell_colors),
  label_x=c(
     2.8,  0.5, -3.9,  0.8, -5.2, -5.5, -5.1,  3.0,
     6.4,  8.0,  8.3,  4.0,  5.5,
    -2.1,  4.9,
   -12.1, -8.5,-13.5, -8.3,-10.4,-13.0, -7.1, -7.4,-12.2
  ),
  label_y=c(
    -3.7, -3.1, -1.9,  3.0,  0.3,  2.4,  4.4,  1.7,
    -1.8, -3.3, -5.3, -7.7,  0.0,
     8.1,  9.1,
    -7.5, -9.3, -5.7, -3.7, -3.0, -3.9, -5.3, -7.1,  9.6
  )
)
label_pos <- label_pos[label_pos$Celltype %in% df$Celltype, ]
label_pos$anchor_x <- label_pos$anchor_y <- NA_real_

# A leader points to a real cell near its group's coordinate median.
for (i in seq_len(nrow(label_pos))) {
  z <- df[df$Celltype == label_pos$Celltype[i], ]
  j <- which.min((z$UMAP_1-median(z$UMAP_1))^2 +
                 (z$UMAP_2-median(z$UMAP_2))^2)
  label_pos$anchor_x[i] <- z$UMAP_1[j]
  label_pos$anchor_y[i] <- z$UMAP_2[j]
}

p <- ggplot(df, aes(UMAP_1, UMAP_2, colour=Celltype)) +
  geom_point(size=1.25, alpha=1, stroke=0) +
  scale_colour_manual(values=cell_colors, breaks=names(cell_colors),
                      drop=FALSE, name=NULL) +
  geom_segment(
    data=label_pos, inherit.aes=FALSE,
    aes(x=anchor_x, y=anchor_y, xend=label_x, yend=label_y),
    colour="#555555", linewidth=0.4
  ) +
  geom_label(
    data=label_pos, inherit.aes=FALSE,
    aes(label_x, label_y, label=Celltype),
    family="Arial", size=4.9, colour="black",
    fill=alpha("white", 0.93), label.size=0,
    label.padding=grid::unit(0.10, "lines")
  ) +
  annotate("segment", x=-14.1, xend=-11.0, y=-10.6, yend=-10.6,
           linewidth=0.8, colour="black",
           arrow=grid::arrow(length=grid::unit(2.8, "mm"), type="closed")) +
  annotate("segment", x=-14.1, xend=-14.1, y=-10.6, yend=-7.7,
           linewidth=0.8, colour="black",
           arrow=grid::arrow(length=grid::unit(2.8, "mm"), type="closed")) +
  annotate("text", x=-12.55, y=-11.25, label="UMAP 1",
           family="Arial", size=6) +
  annotate("text", x=-14.7, y=-9.15, label="UMAP 2", angle=90,
           family="Arial", size=6) +
  coord_equal(xlim=c(-15.5, 10), ylim=c(-11.7, 11.2), clip="off") +
  theme_void(base_family="Arial", base_size=20) +
  theme(
    legend.position="left",
    legend.text=element_text(size=16, colour="black"),
    legend.key.height=grid::unit(7.3, "mm"),
    legend.key.width=grid::unit(7, "mm"),
    legend.box.margin=margin(0, 12, 0, 0),
    plot.margin=margin(10, 14, 10, 10),
    plot.background=element_rect(fill="white", colour=NA)
  ) +
  guides(colour=guide_legend(
    ncol=1, override.aes=list(size=3.5, alpha=1)
  ))

p
