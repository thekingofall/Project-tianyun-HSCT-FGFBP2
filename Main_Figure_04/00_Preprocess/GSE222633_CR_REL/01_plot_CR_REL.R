library(ggplot2)
if (length(list.files("figures/png")) || length(list.files("figures/pdf")) ||
    file.exists("data/plot_ready/plotted_patient_values.csv"))
  stop("Existing outputs retained; use a new version directory")

# Run from the delivered CR_REL_Teff version folder.
d <- read.csv("data/plot_ready/sample_Teff_expression.csv",
              stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(d[c("sample", "gene")]),
          all(d$outcome == ifelse(d$Pair == "D", "REL", "CR")))
d$timepoint <- factor(d$timepoint,
                     levels = c("Pre_GCSF", "Post_GCSF", "D90", "D180"))
d$x <- as.integer(d$timepoint)
d$Pair <- factor(d$Pair, levels = LETTERS[1:5])
d$outcome <- factor(d$outcome, levels = c("CR", "REL"))
# One common gene-wise scale across all 17 sample means, before grouping.
d$z_score <- ave(d$mean_log_expression, d$gene,
                 FUN = function(v) (v - mean(v)) / sd(v))
# Deterministic horizontal offsets separate overlapping patient points.
offsets <- c(A=-0.10, B=-0.035, C=0.035, D=0, E=0.10)
d$x_point <- d$x + offsets[as.character(d$Pair)]
pal <- read.delim("palette_registry.tsv", stringsAsFactors=FALSE)
colors <- setNames(pal$hex, pal$outcome)
summary_data <- aggregate(z_score ~ gene + outcome + x, d, mean)
summary_data$n_pairs <- vapply(seq_len(nrow(summary_data)), function(i) {
  with(summary_data[i, ], sum(d$gene == gene & d$outcome == outcome & d$x == x))
}, integer(1))
write.csv(d, "data/plot_ready/plotted_patient_values.csv", row.names=FALSE)
write.csv(summary_data, "data/plot_ready/plotted_group_means.csv", row.names=FALSE)

for (g in c("FGFBP2", "ADGRG1")) {
  points <- subset(d, gene == g)
  means <- subset(summary_data, gene == g)
  p <- ggplot() +
    geom_line(data=means,
              aes(x=x, y=z_score, color=outcome, group=outcome),
              linewidth=1.15) +
    geom_point(data=points,
               aes(x=x_point, y=z_score, color=outcome, shape=Pair),
               size=3.3, stroke=0.8) +
    scale_color_manual(values=colors, drop=FALSE) +
    scale_shape_manual(values=c(A=16, B=15, C=18, D=17, E=8),
                       labels=paste("Pair", LETTERS[1:5]), drop=FALSE) +
    scale_x_continuous(breaks=1:4,
                       labels=c("Pre-GCSF", "Post-GCSF", "D90", "D180"),
                       limits=c(0.75,4.25)) +
    scale_y_continuous(breaks=-2:2, limits=c(-2.1,2.1)) +
    labs(x=NULL, y=paste0(g, " expression (z-score)"),
         color=NULL, shape=NULL) +
    theme_classic(base_family="Arial", base_size=17) +
    theme(axis.title.y=element_text(face="bold", margin=margin(r=12)),
          axis.text=element_text(color="black"),
          axis.text.x=element_text(angle=35, hjust=1),
          axis.line=element_line(linewidth=0.85, color="black"),
          axis.ticks=element_line(linewidth=0.7),
          legend.position="right", legend.box="vertical",
          legend.text=element_text(size=14),
          plot.margin=margin(12,14,10,10)) +
    guides(color=guide_legend(order=1, override.aes=list(linewidth=1.2)),
           shape=guide_legend(order=2, override.aes=list(color="#4D4D4D")))
  png_file <- paste0("figures/png/",g,"_Teff_CR_REL_zscore.png")
  pdf_file <- paste0("figures/pdf/",g,"_Teff_CR_REL_zscore.pdf")
  if (file.exists(png_file) || file.exists(pdf_file))
    stop("Existing figure retained; run in a new version directory")
  print(p)
  ggsave(png_file,p,width=7,height=4.8,dpi=300,bg="white",type="cairo")
  ggsave(pdf_file,p,width=7,height=4.8,device=cairo_pdf,family="Arial")
}
writeLines(capture.output(sessionInfo()), "logs/sessionInfo.txt")
cat("DONE: two Teff CR/REL figures\n")
