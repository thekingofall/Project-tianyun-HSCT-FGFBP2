library(ggplot2)
library(grid)
# Run from this delivered version folder. Each rerun creates a new folder.
out <- Sys.getenv("FIGURE_OUTPUT_ROOT", unset=paste0("visualization_",format(Sys.time(),"%Y%m%d_%H%M%S")))
for (sub in c("figures/png","figures/pdf","data/plot_ready","logs"))
  dir.create(file.path(out,sub), recursive=TRUE, showWarnings=FALSE)
pngfile <- file.path(out,"figures/png/FGFBP2_Teff_CR_REL_trajectories_and_means.png")
pdffile <- file.path(out,"figures/pdf/FGFBP2_Teff_CR_REL_trajectories_and_means.pdf")
if (file.exists(pngfile) || file.exists(pdffile)) stop("Existing figures retained; use a new output directory")
d <- read.csv("data/plot_ready/FGFBP2_sample_expression.csv",stringsAsFactors=FALSE)
stopifnot(nrow(d)==17L, !anyDuplicated(d$sample), all(d$outcome==ifelse(d$Pair=="D","REL","CR")))
d$x <- match(d$timepoint,c("Pre_GCSF","Post_GCSF","D90","D180"))
d$z <- as.numeric(scale(d$mean_log_expression))
d$Pair <- factor(d$Pair,levels=LETTERS[1:5])
d$outcome <- factor(d$outcome,levels=c("CR","REL"))
d <- d[order(d$Pair,d$x),]
registry <- read.delim("palette_registry.tsv")
cols <- setNames(registry$hex,registry$label)
shapes <- c(A=16,B=15,C=18,D=17,E=8)
# Connect adjacent visits only. Post-GCSF to D90 links donor and recipient.
segs <- do.call(rbind,lapply(split(d,d$Pair),function(a) {
 if(nrow(a)<2) return(NULL)
 b <- a[-nrow(a),]; b$xend <- a$x[-1]; b$yend <- a$z[-1]
 b <- b[b$xend-b$x==1,]
 b$bridge <- ifelse(b$x==2,"Donor to recipient","Same participant")
 b
}))
m <- aggregate(mean_log_expression~outcome+x,d,mean)
m$n <- vapply(seq_len(nrow(m)),function(i)sum(d$outcome==m$outcome[i]&d$x==m$x[i]),integer(1))
d$xd <- d$x + ifelse(d$outcome=="CR",-.15,.15)
d$xp <- d$xd + c(A=-.05,B=-.016,C=.016,D=0,E=.05)[as.character(d$Pair)]
m$xd <- m$x + ifelse(m$outcome=="CR",-.15,.15)
write.csv(d,file.path(out,"data/plot_ready/plotted_samples.csv"),row.names=FALSE)
write.csv(m,file.path(out,"data/plot_ready/stage_means.csv"),row.names=FALSE)
base <- theme_classic(base_family="Arial",base_size=16)+
 theme(axis.text=element_text(colour="black"),axis.title.y=element_text(face="bold",margin=margin(r=10)),
       axis.line=element_line(colour="black",linewidth=.8),axis.ticks=element_line(linewidth=.7),
       legend.position="right",legend.text=element_text(size=14),
       plot.margin=margin(12,8,8,8))
xt <- scale_x_continuous(breaks=1:4,limits=c(.6,4.4),expand=c(0,0))
p1 <- ggplot(d,aes(x=x,y=z,colour=outcome))+
 geom_segment(data=segs,aes(xend=xend,yend=yend,linetype=bridge),linewidth=.8,show.legend=FALSE)+
 geom_point(aes(shape=Pair),size=3.7,stroke=.85)+
 scale_colour_manual(values=cols)+
 scale_shape_manual(values=shapes,labels=paste("Pair",LETTERS[1:5]),drop=FALSE)+
 scale_linetype_manual(values=c("Same participant"="solid","Donor to recipient"="dashed"))+
 xt+scale_y_continuous(breaks=-1:2,limits=c(-1.5,2.15))+
 annotate("text",x=c(1.5,3.5),y=2.05,label=c("Donor","Recipient"),family="Arial",size=4.8)+
 labs(x=NULL,y="FGFBP2 expression (z-score)",colour=NULL,shape=NULL)+base+
 theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
 guides(colour=guide_legend(order=1,override.aes=list(shape=16,size=4)),
        shape=guide_legend(order=2,override.aes=list(colour="#444444",size=3.7)))
p2 <- ggplot(d,aes(x=xp,y=mean_log_expression,colour=outcome))+
 geom_point(aes(shape=Pair),size=3.7,stroke=.85)+
 geom_segment(data=m,aes(x=xd-.085,xend=xd+.085,y=mean_log_expression,yend=mean_log_expression),
              linewidth=1.6,inherit.aes=FALSE,colour=cols[as.character(m$outcome)])+
 scale_colour_manual(values=cols)+scale_shape_manual(values=shapes)+
 scale_x_continuous(breaks=1:4,limits=c(.6,4.4),expand=c(0,0),
  labels=c("Pre-GCSF\n4 CR / 0 REL","Post-GCSF\n4 CR / 1 REL","D90\n3 CR / 1 REL","D180\n3 CR / 1 REL"))+
 scale_y_continuous(limits=c(0,2.6),breaks=c(0,1,2),expand=c(0,0))+
 labs(x=NULL,y="Mean log-normalized expression")+base+
 theme(legend.position="none",axis.text.x=element_text(size=12,lineheight=1.15))
g <- ggplotGrob(p1)
legend <- g$grobs[[which(vapply(g$grobs,function(x)x$name,character(1))=="guide-box")]]
g1 <- ggplotGrob(p1+theme(legend.position="none"));g2 <- ggplotGrob(p2)
w <- unit.pmax(g1$widths,g2$widths);g1$widths <- w;g2$widths <- w
# A shared grid guarantees the four stage centers align across panels.
draw <- function() {
 grid.newpage();pushViewport(viewport(layout=grid.layout(2,2,widths=unit(c(1,.20),"null"),heights=unit(c(1.12,1),"null"))))
 pushViewport(viewport(layout.pos.row=1,layout.pos.col=1));grid.draw(g1);popViewport()
 pushViewport(viewport(layout.pos.row=2,layout.pos.col=1));grid.draw(g2);popViewport()
 pushViewport(viewport(layout.pos.row=1:2,layout.pos.col=2));grid.draw(legend);popViewport(2)
}
png(pngfile,width=8.5,height=7.2,units="in",res=300,type="cairo",bg="white");draw();dev.off()
cairo_pdf(pdffile,width=8.5,height=7.2,family="Arial");draw();dev.off()
writeLines(capture.output(sessionInfo()),file.path(out,"logs/sessionInfo.txt"))
cat("DONE FGFBP2; 17 samples; 5 Pairs; figures exported\n")
