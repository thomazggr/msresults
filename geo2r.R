# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
################################################################
#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)
#library(ggplot2)

Sys.setenv(VROOM_CONNECTION_SIZE = "500000")

# load series and platform data from GEO

run.geo <- function(gse, annot_gpl = FALSE, gpl, samples){
  gset <- getGEO(gse, GSEMatrix =TRUE, AnnotGPL=annot_gpl)
  if (length(gset) > 1) idx <- grep(gpl, attr(gset, "names")) else idx <- 1
  gset <- gset[[idx]]

  # make proper column names to match toptable 
  fvarLabels(gset) <- make.names(fvarLabels(gset))

  # group membership for all samples
  gsms <- samples
  sml <- strsplit(gsms, split="")[[1]]

  # filter out excluded samples (marked as "X")
  sel <- which(sml != "X")
  sml <- sml[sel]
  gset <- gset[ ,sel]

  # log2 transformation
  ex <- exprs(gset)
  qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
            (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { ex[which(ex <= 0)] <- NaN
    exprs(gset) <- log2(ex) }

  # assign samples to groups and set up design matrix
  gs <- factor(sml)
  groups <- make.names(c("CT","NASH"))
  levels(gs) <- groups
  gset$group <- gs
  design <- model.matrix(~group + 0, gset)
  colnames(design) <- levels(gs)

  fit <- lmFit(gset, design)  # fit linear model

  # set up contrasts of interest and recalculate model coefficients
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)

  # compute statistics and table of top significant genes
  fit2 <- eBayes(fit2, 0.01)
  tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

  # # tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
  write.table(tT, file=paste0("~/Documents/Projects/msresults/", gse, ".tsv"), row.names=F, sep="\t")

  # # Visualize and quality control test results.
  # # Build histogram of P-values for all genes. Normal test
  # # assumption is that most genes are not differentially expressed.
  # tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
  # png(file=paste0(gse, "_histogram.png"), width=840, height=720)
  # hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
    # ylab = "Number of genes", main = "P-adj value distribution")
  # dev.off()
  
  # # # summarize test results as "up", "down" or "not expressed"
  # dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05)

  # # # Venn diagram of results
  # # vennDiagram(dT, circle.col=palette())

  # # # create Q-Q plot for t-statistic
  # # t.good <- which(!is.na(fit2$F)) # filter out bad probes
  # # qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

  # # # volcano plot (log P-value vs log fold change)
  # colnames(fit2) # list contrast names
  # ct <- 1        # choose contrast of interest
  
  # png(file=paste0(gse, "_volcano.png"), width=720, height=840)
  # limma::volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
    # highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))
  # dev.off()

  # # MD plot (log fold change vs mean log expression)
  # # highlight statistically significant (p-adj < 0.05) probes
  # png(file=paste0(gse, "_mdplot.png"), width=840, height=720)
  # limma::plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
  # abline(h=0)
  # dev.off()

  # # ################################################################
  # # # General expression data analysis
  # ex <- exprs(gset)

  # # # box-and-whisker plot
  # ord <- order(gs)  # order samples by group
  # palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
            # "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
  # par(mar=c(7,4,2,1))
  # title <- paste (gse, "/", annotation(gset), sep ="")
  # png(file=paste0(gse, "_box_and_whisker.png"), width=840, height=720)
  # boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
  # legend("topleft", groups, fill=palette(), bty="n")
  # dev.off()

  # # # expression value distribution
  # par(mar=c(4,4,2,1))
  # title <- paste ("GSE33814", "/", annotation(gset), " value distribution", sep ="")
  # png(file=paste0(gse, "_expression_dist.png"), width=840, height=720)
  # limma::plotDensities(ex, group=gs, main=title, legend ="topright")
  # dev.off()

  # # UMAP plot (dimensionality reduction)
  # ex <- na.omit(ex) # eliminate rows with NAs
  # ex <- ex[!duplicated(ex), ]  # remove duplicates
  # ump <- umap(t(ex), n_neighbors = 11, random_state = 123)
  # par(mar=c(3,3,2,6), xpd=TRUE)
  # plot(ump$layout, main="UMAP plot, nbrs=11", xlab="", ylab="", col=gs, pch=20, cex=1.5)
  # legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
  # col=1:nlevels(gs), title="Group", pt.cex=1.5)
  # library("maptools")  # point labels without overlaps
  # pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)

  # # mean-variance trend, helps to see if precision weights are needed
  # plotSA(fit2, main="Mean variance trend, GSE33814")
}

datasets <- list(
  c("GSE33814", "1XXX1X111111111XX1XXXXXX000XX0XXX00000X00X00", TRUE, "GPL6884"),
  c("GSE89632", "X1XXXXXX11X11XX111X1X111XXXX1X11111X10000000000000000X000X00000", FALSE, "GPL14951"),
  c("GSE37031", "111111110000000", FALSE, "GPL14877"),
  c("GSE48452", "000001X0001XXXXXXXXXXXXXX11X0XX11XXX1X1X0XXXX1X10XXXXX0X1XXXXXX1X11XXXXXX", TRUE, "GPL11532"),
  c("GSE63067", "XX1111111110000000", TRUE, "GPL570")
)

df_datasets <- as.data.frame(do.call(rbind, datasets))

colnames(df_datasets) <- c("gse", "samples", "annot_gpl", "gpl")

for(i in 1:nrow(df_datasets)) {
  row <- df_datasets[i, ]

  run.geo(gse=row$gse, annot_gpl=row$annot_gpl, gpl=row$gpl, samples=row$samples)
}
