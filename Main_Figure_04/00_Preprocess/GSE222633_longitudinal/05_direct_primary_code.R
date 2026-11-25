# Set working directory to the delivered version folder.
# These are the exact compact data used for the delivered primary figures.
library(data.table)
library(ggplot2)
d <- fread("data/plot_ready/Teff_gene_expression.csv")
pal <- c(A="#3B6C80", B="#B47A4E", C="#68875A", D="#A34F65", E="#7C6A98")
d[, x := match(timepoint, c("Pre_GCSF", "Post_GCSF", "D90", "D180"))]
d[, Pair := factor(Pair, levels=LETTERS[1:5])]
for (g in c("FGFBP2", "ADGRG1")) {
  # Full grids prevent lines bridging an unobserved timepoint.
  q <- merge(CJ(Pair=factor(LETTERS[1:5], levels=LETTERS[1:5]), x=1:4),
             d[gene==g, .(Pair, x, mean_log_expression)], by=c("Pair", "x"), all.x=TRUE)
  setorder(q, Pair, x)
  q[, `:=`(xend=shift(x, type="lead"), yend=shift(mean_log_expression, type="lead")), by=Pair]
  seg <- q[is.finite(mean_log_expression) & is.finite(yend)]
  seg[, bridge := ifelse(x==2, "Donor to recipient", "Within person")]
  p <- ggplot(q, aes(x, mean_log_expression, color=Pair)) +
    geom_segment(data=seg, aes(xend=xend, yend=yend, linetype=bridge), linewidth=.85,
                 show.legend=c(linetype=FALSE)) +
    geom_segment(data=seg[Pair=="D"], aes(xend=xend, yend=yend, linetype=bridge),
                 linewidth=1.35, show.legend=FALSE) +
    geom_point(aes(shape=Pair), size=3.3, na.rm=TRUE) +
    scale_color_manual(values=pal, labels=paste("Pair", LETTERS[1:5])) +
    scale_shape_manual(values=c(16,17,15,18,8), labels=paste("Pair", LETTERS[1:5])) +
    scale_linetype_manual(values=c("Donor to recipient"="dashed", "Within person"="solid")) +
    scale_x_continuous(breaks=1:4, labels=c("Pre","Post","D90","D180"), limits=c(.85,4.15)) +
    labs(x=NULL, y=paste0(g," expression\n(mean log-normalized RNA)"), color=NULL, shape=NULL) +
    theme_classic(base_family="Arial", base_size=17) +
    theme(axis.title=element_text(face="bold"), axis.text=element_text(color="black"),
          axis.line=element_line(linewidth=.85), legend.position="right")
  print(p)
  out <- paste0(g, "_Teff_trajectory_redraw")
  stopifnot(!file.exists(paste0(out,".png")), !file.exists(paste0(out,".pdf")))
  ggsave(paste0(out,".png"), p, width=8, height=5, dpi=300, bg="white")
  ggsave(paste0(out,".pdf"), p, width=8, height=5, device=cairo_pdf, family="Arial")
}
