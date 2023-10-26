#' ---
#' title: "Dimensionality reduction using zinbwave, initially designed for AD only, to get a clean GammaDelta cluster"
#' author: "gabrielet"
#' output: html_document
#' date: '`r format(Sys.Date(), "%d %B, %Y")`'
#' ---
#' 
#' ```{r setup, include=FALSE}
#' knitr::opts_chunk$set(
#'   tidy=TRUE,
#'   tidy.opts=list(width.cutoff=120),
#'   message=FALSE,
#'   warning=FALSE,
#'   eval=TRUE
#' )
#' ```

# Load libraries
library("bluster")
library("cluster")
library("compositions")
library("cowplot")
library("ggplot2")
library("gridExtra")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library("BiocParallel")
library("zinbwave")
library("clusterExperiment")

# set paths
rawPath <- "/home/gabriele/work/cbmc/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"

saveIn <- paste0(rawPath, "ANALYSIS_FINAL_441/")
saveImg <- paste0(saveIn, "figures/")

print(paste0("analysing exp ", expN))

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# load clean QCed data to export only
cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

# get best 1000 genes
perC <- 1000
topGMat <- as.matrix(assays(cleanSce)$counts)
topGMat <- makeFilterStats(topGMat, filterStats="var", transFun=log1p)
topGMat <- filterData(topGMat, percentile=perC, filterStats="var")

sceTop1000 <- cleanSce[rownames(topGMat),]
topGMat$condition <- sceTop1000$condition
topGMat$pruned_main <- sceTop1000$pruned_main
topGMat$pruned_fine <- sceTop1000$pruned_fine
assayNames(topGMat) <- "counts"

# export topGMat to save some ram, then restart R, re-load libraries and paths then 

saveRDS(object=topGMat, file=paste0(saveIn, "topGMat_data_", expN, ".Rds"))

################################################################ RE-START FROM HERE

####################get best 1000 genes
###################perC <- 1000

###################topGMat <- readRDS(paste0(saveIn, paste0("topGMat_data_", expN, ".Rds")))

#and finally run zinbwave

# best parameter combination is perC = 1000 and kZinb = 20

# it may happen that a warning arise at the end of the computation stating that:
# Warning message: In simpute.als(x, J, thresh, lambda, maxit, trace.it, warm.start, : Convergence not achieved by 100 iterations
# this is not a problem, according to this page:
# https://github.com/drisso/zinbwave/issues/9

kZinb <- 20

########################################set.seed(131)

# run zinbwave
########################################print(system.time(zbWaved <- zinbwave(Y=topGMat, K=kZinb, X="~condition", which_assay="counts", residuals=TRUE, normalizedValues=TRUE, verbose=TRUE, BPPARAM=BiocParallel::SerialParam())))

# export the just zinbwaved object
########################################saveRDS(object=zbWaved, file=paste0(saveIn, "just_zbwaved_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))

######################################### naming assays
########################################assays(zbWaved) <- assays(zbWaved)[c("normalizedValues", "counts", "residuals")]

######################################### using normalised values as logcounts which will be used by findMarkers, later
########################################assays(zbWaved)$logcounts <- assays(zbWaved)$normalizedValues

# also, consider that using normalised values as logcounts will make sure these value will be used by findMarkers









# LOAD OBJECT
zbWaved <- readRDS(paste0(saveIn, "full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))

# SET SEED
set.seed(131)

# then run tsne
zbWaved <- runTSNE(zbWaved, dimred="zinbwave", perplexity=30)

# and plot
png(paste0(saveImg, "02_01_zinbwave_tSNE_", expN, "_K_", kZinb, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	scater::plotTSNE(zbWaved, colour_by="condition")
dev.off()

# and plot
png(paste0(saveImg, "02_02_zinbwave_tSNE_condition_", expN, "_K_", kZinb, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
grid.arrange(
	scater::plotReducedDim(zbWaved, dimred="TSNE", colour_by="condition") + theme(legend.position="right") + labs(title="By condition") + guides(color=guide_legend("Condition")),
	scater::plotReducedDim(zbWaved, dimred="TSNE", colour_by="pruned_main") + theme(legend.position="right") + labs(title="By cell type") + guides(color=guide_legend("Cell type")) ,
	ncol=2)
dev.off()

# and export a summary table
summ <- table(zbWaved$pruned_main, zbWaved$condition)
summarY <- data.frame(rownames(summ), summ[, 1], summ[, 2])
colnames(summarY) <- c("cType", "3xTG", "WT")
write.table(summarY, paste0(saveImg, "cells_X_condition_", expN, "_K_", kZinb, ".csv"), quote=F, sep="\t", row.names=F)

# finally, export the data
saveRDS(object=zbWaved, file=paste0(saveIn, "full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))
