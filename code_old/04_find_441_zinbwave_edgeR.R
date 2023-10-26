#' ---
#' title: "Analysing GammaDelta cells in AD exps"
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
library("DropletUtils")
library("ggplot2")
library("gridExtra")
library("pheatmap")
library("plyr")
library("limma")
library("clusterProfiler")
library("org.Mm.eg.db")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library("openxlsx")
library("fgsea")
library("pathview")
library("ReactomePA")
library("DESeq2")
library("enrichplot")
library("DOSE")
library("HelpersMG")
library("gg.gap")
library("ddpcr")
library("data.table")
library("ggpubr")
library("edgeR")

# import plotting functions
source("00_functions.R")

# initialise label
lbl <- 1

# set path
rawPath <- "/home/gabriele/Desktop/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"

cellType <- c("Neutrophils", "CD4", "Tgd")
cellLab <-  c("Neutrophils", "CD4", "Tgd")

# analyse TDG only, for now
cT <- 3

print("")
print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
print("")

# then do everything
saveIn <- paste0(rawPath, "results_ribosomal_441_no_ribosomes_2_luglio_sera_senza_RIKEN/")
saveUpDown <- paste0(saveIn, "enrichment_", cellLab[cT], "/UpDown/")
saveImg <- paste0(saveIn, "figures/analysis_", cellLab[cT], "/")

# check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveUpDown), TRUE, dir.create(saveUpDown, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

# NOTE: we are using QCed data file which logcounts were normalised using sizeFactors
# but this data file contains also the information from the SLM clustering
# performed using the zinbwaved data, in order to have a good clustering but also all
# the genes, instead of the best 1000 used for zinbwave
perC <- 1000
kZinb <- 20
clusteredSce <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, "_withRibosomes.Rds"))

# create fine_to_main labels as it was done for the zbWaved object in 03
clusteredSce$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(clusteredSce$pruned_fine, "\\("), `[[`, 1)))

# plotting heatmap with cell type
tab_fine_to_main <- table(Assigned=clusteredSce$fine_to_main, Cluster=clusteredSce$slm_0.5)

# get the table ready to be exported
eTab <- as.data.frame.matrix(tab_fine_to_main)
eTabRN <- cbind.data.frame(cellType=rownames(eTab), eTab)
eTabCN <- rbind.data.frame(c("cluster", colnames(eTab)), eTabRN)

# export the table
write.table(eTabCN, paste0(saveImg, "cluster_content_exp_", expN, "_", cellLab[cT], ".csv"), quote=F, row.names=F, col.names=F, sep="\t")
write.xlsx(eTabCN, paste0(saveImg, "cluster_content_exp_", expN, "_", cellLab[cT], ".xlsx"), colNames=F)

tab_full_labels <- table(Assigned=clusteredSce$pruned_fine, Cluster=clusteredSce$slm_0.5)

# get the full table ready to be exported
etabF <- as.data.frame.matrix(tab_full_labels)
eTabrnF <- cbind.data.frame(cellType=rownames(etabF), etabF)
eTabcnF <- rbind.data.frame(c("cluster", colnames(etabF)), eTabrnF)

# export the full table
write.table(eTabcnF, paste0(saveImg, "cluster_content_full_labels_exp_", expN, "_", cellLab[cT], ".csv"), quote=F, row.names=F, col.names=F, sep="\t")
write.xlsx(eTabcnF, paste0(saveImg, "cluster_content_full_labels_exp_", expN, "_", cellLab[cT], ".xlsx"), colNames=F)

# get data set conditions, i.e. WT and 3xTG
conD <- unique(clusteredSce$condition)

# search for cellType
cTPos <- grep(cellType[cT], clusteredSce$pruned_fine)

print(paste0(cellLab[cT], " are found in ", length(unique(clusteredSce$slm_0.5[cTPos])), " different clusters"))

fR <- table(clusteredSce$slm_0.5[cTPos])
clustS <- names(fR[which(fR!=0)])

print(fR[which(fR!=0)])

# Cluster composition plots using main labels

df <- data.frame(table(clusteredSce$slm_0.5, clusteredSce$fine_to_main, clusteredSce$condition))

colnames(df) <- c("Cluster", "CellType", "Condition", "Frequency")
lims <- levels(factor(clusteredSce$fine_to_main))
lims <- lims[order(table(clusteredSce$fine_to_main))]

# plot!
p <- plot_membership(df, cond=conD[2], col=colorsPal[2])
# update label
lbl <- lbl + 1
png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(p)
dev.off()

p <- plot_membership(df, cond=conD[1], col=colorsPal[6])
# update label
lbl <- lbl + 1
png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	plot(p)
dev.off()

print(paste0("######## SEARCHING FOR ", cellType[cT], " IN THE WHOLE DATASET ########"))

# search for cellType
cTPos <- grep(cellType[cT], clusteredSce$pruned_fine)

# and subset the original single cell experiment using the name of the cellType cells.
cTSubSet <- clusteredSce[, cTPos]

# set WT as a reference, to compute fold change
cTSubSet$condition <- relevel(as.factor(cTSubSet$condition), ref="WT")

# set test
tst <- "findM"

if (tst == "deseq") {

	cat("using deseq\n")
	
	design <- stats::model.matrix(object = ~ condition, data = colData(cTSubSet))

	dds <- DESeqDataSetFromMatrix(countData = counts(cTSubSet),
		                      colData = colData(cTSubSet),
		                      design = design)
	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	dds <- computeSumFactors(dds, min.mean=0.5)
	dds <- logNormCounts(dds)

	# get the just computed sizeFactors
	normFacts <- dds$sizeFactor

	# renormalize: multiply to 1
	normFacts <- normFacts/exp(mean(log(normFacts)))

	sizeFactors(dds) <- normFacts

	library("glmGamPoi")

	dds <- DESeq(dds, useT = TRUE, minmu = 1e-6, minReplicatesForReplace = Inf, sfType = "ratio", fitType = "glmGamPoi", test = "LRT", full = model.matrix(~ 1 + condition, data = colData(cTSubSet)), reduced = model.matrix(~ 1, colData(cTSubSet)))

	res <- results(dds)

	statInfo <- as.data.frame(res)

	# get significant
	siG <- statInfo[statInfo$padj <= 0.05, ]
	dim(siG[siG$log2FoldChange >= 0.5, ])
	dim(siG[siG$log2FoldChange < (-0.5), ])
	dim(siG[siG$log2FoldChange >= 0, ])
	dim(siG[siG$log2FoldChange < 0, ])

} else if (tst == "edger") {

	cat("using edgeR\n")

	# remove genes with too many zeros
#	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05,]
	
	# nuovo filtro
	cTSubSet <- cTSubSet[rowSums(counts(cTSubSet) > 0) > 10,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)

	# get the just computed sizeFactors
	normFacts <- cTSubSet$sizeFactor

	# renormalize: multiply to 1
	normFacts <- normFacts/exp(mean(log(normFacts)))

	# now run edgeR WITH glmQLFTest

	# create DGE object
	y <- DGEList(counts = assays(cTSubSet)$counts, group = cTSubSet$condition, norm.factors = normFacts)

	# set parameters for computing glmQLFTest
	design <- stats::model.matrix(object = ~ condition, data = colData(cTSubSet))
	y <- edgeR::estimateDisp(y = y, design = design)
	dispEsts <- edgeR::getDispersion(y)
	glmFit <- edgeR::glmQLFit(y = y, dispersion = dispEsts, robust = FALSE, design = design)
	glmRes <- edgeR::glmQLFTest(glmFit, coef = 2)
	statInfo <- glmRes[["table"]]
	pval <- statInfo[, "PValue"]
	pValMat <- data.frame("rawP" = pval, "adjP" = stats::p.adjust(pval, "BH"))
	rownames(pValMat) <- rownames(statInfo)

	# finalise tables
	allRes <- cbind.data.frame(gName=rownames(statInfo), statInfo, adjP=pValMat$adjP)
	finalUns <- allRes[which(allRes$adjP <= 0.05), ]

	# sort by logfoldchange
	sortedFinal <- finalUns[order(finalUns$logFC, decreasing=T), ]

	# and export
	write.table(sortedFinal, paste0(saveUpDown, "up_down_EDGER_", cellLab[cT], "_CELLS_", expN, ".csv"), col.names=T, row.names=F, sep="\t")
	write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_EDGER_", cellLab[cT], "_CELLS_", expN, ".xlsx"))

} else if (tst == "findM") {

	cat("using findMarkers\n")

	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05,]
	
	# nuovo filtro
#	cTSubSet <- cTSubSet[rowSums(counts(cTSubSet) > 0) > 10,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)

	cTSce <- list()
	conD <- unique(clusteredSce$condition)
	# then get WT cellType
	cTSce[[conD[1]]] <- cTSubSet[, which(cTSubSet$condition==conD[1])]
	# and 3xTG
	cTSce[[conD[2]]] <- cTSubSet[, which(cTSubSet$condition==conD[2])]


	markersG <- findMarkers(cTSubSet, groups=cTSubSet$condition, pval.type="any", direction="any")

	cC <- "3xTG"
	chosen <- markersG[[cC]]
	selfMade <- chosen[chosen$FDR < 0.05, ]
	# and export data
	finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)
	
	# sort by logfoldchange
	sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

	# and export
	write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_", cellLab[cT], "_CELLS_", expN, ".csv"), col.names=T, row.names=F, sep="\t")
	write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_", cellLab[cT], "_CELLS_", expN, ".xlsx"))
}

####################### print some violins

# FOR DOWN REG

saveExpression <- paste0(saveIn, "violins_DOWNReg/")
ifelse(dir.exists(saveExpression), TRUE, dir.create(saveExpression, recursive=T))

genesToPlot <- sortedFinal$gName[sortedFinal$logFC < 0]

# loop through the genes of interest
for (gTP in genesToPlot) {
	if (gTP %in% rownames(cTSubSet)) {
		if (all(logcounts(cTSubSet[rownames(cTSubSet) %in% gTP, ]) == 0)) {
			print(paste0("all logcounts equal to zero. skipping ", gTP))
		} else {
			# plot the logcounts for cellType cells only, as violins, to compare the two conditions
			lbl <- lbl + 1
			p <- plotExpression(cTSubSet, x="condition", colour_by="condition", features=gTP, show_median=T) + stat_compare_means(method = "wilcox.test", label.x=1.5)
			png(paste0(saveExpression, "04_0", lbl, "_comparing_", gTP, "_by_condition_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
				plot(p)
			dev.off()
		}
	} else {
		print(paste0("gene ", gTP, " not found!"))
	}
}

# AND UP REG

saveExpression <- paste0(saveIn, "violins_UPReg/")
ifelse(dir.exists(saveExpression), TRUE, dir.create(saveExpression, recursive=T))

genesToPlot <- sortedFinal$gName[sortedFinal$logFC >= 0]

# loop through the genes of interest
for (gTP in genesToPlot) {
	if (gTP %in% rownames(cTSubSet)) {
		if (all(logcounts(cTSubSet[rownames(cTSubSet) %in% gTP, ]) == 0)) {
			print(paste0("all logcounts equal to zero. skipping ", gTP))
		} else {
			# plot the logcounts for cellType cells only, as violins, to compare the two conditions
			lbl <- lbl + 1
			p <- plotExpression(cTSubSet, x="condition", colour_by="condition", features=gTP, show_median=T) + stat_compare_means(method = "wilcox.test", label.x=1.5)
			png(paste0(saveExpression, "04_0", lbl, "_comparing_", gTP, "_by_condition_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
				plot(p)
			dev.off()
		}
	} else {
		print(paste0("gene ", gTP, " not found!"))
	}
}
