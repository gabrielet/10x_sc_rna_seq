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
rawPath <- "/home/gabriele/work/cbmc/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"
 
saveIn <- paste0(rawPath, "ANALYSIS_FINAL_441/")

# NOTE: we are using QCed data file which logcounts were normalised using sizeFactors
# but this data file contains also the information from the SLM clustering
# performed using the zinbwaved data, in order to have a good clustering but also all
# the genes, instead of the best 1000 used for zinbwave
perC <- 1000
kZinb <- 20	

# get normally clustered data
clusteredSce <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))

# create fine_to_main labels as it was done for the zbWaved object in 03
clusteredSce$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(clusteredSce$pruned_fine, "\\("), `[[`, 1)))

# select cell types
cellType <- c("T.8|T.CD8", "Neutrophils", "Tgd")#, "Monocytes", "Macrophages", "Microglia", "B cells", "ILC", "NK cells", "NKT", "Mast")
cellLab <-  c("CD8", "Neutrophils", "Tgd")#, "Monocytes", "Macrophages", "Microglia", "B_cells", "ILC", "NK_cells", "NKT", "Mast")

# loop through cell type
for (cT in seq(1, length(cellLab), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
	print("")
	
	# set paths
	saveUpDown <- paste0(saveIn, "enrichment_", cellLab[cT], "/UpDown/")
	saveImg <- paste0(saveIn, "figures/analysis_", cellLab[cT], "/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveUpDown), TRUE, dir.create(saveUpDown, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))

	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")
	
	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=clusteredSce$fine_to_main, Cluster=clusteredSce$slm_0.5, Condition=clusteredSce$condition)

	# get 3x table ready to be exported
	eTabTX <- as.data.frame.matrix(tab_fine_to_main[, , "3xTG"])
	eTabRNTX <- cbind.data.frame(cellType=rownames(eTabTX), eTabTX)
	eTabCNTX <- rbind.data.frame(c("cluster", colnames(eTabTX)), eTabRNTX)

	# same for WT
	eTabWT <- as.data.frame.matrix(tab_fine_to_main[, , "WT"])
	eTabRNWT <- cbind.data.frame(cellType=rownames(eTabWT), eTabWT)
	eTabCNWT <- rbind.data.frame(c("cluster", colnames(eTabWT)), eTabRNWT)

	# export both tables
	write.xlsx(eTabCNTX, paste0(saveImg, "TX_cluster_content_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabCNWT, paste0(saveImg, "WT_cluster_content_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=F)

	##################### FULL

	tab_full_labels <- table(Assigned=clusteredSce$pruned_fine, Cluster=clusteredSce$slm_0.5, Condition=clusteredSce$condition)

	# get 3xTG full table ready to be exported
	etabFTX <- as.data.frame.matrix(tab_full_labels[, , "3xTG"])
	eTabrnFTX <- cbind.data.frame(cellType=rownames(etabFTX), etabFTX)
	eTabcnFTX <- rbind.data.frame(c("cluster", colnames(etabFTX)), eTabrnFTX)

	# and WT
	etabFWT <- as.data.frame.matrix(tab_full_labels[, , "WT"])
	eTabrnFWT <- cbind.data.frame(cellType=rownames(etabFWT), etabFWT)
	eTabcnFWT <- rbind.data.frame(c("cluster", colnames(etabFWT)), eTabrnFWT)

	# export the full table
	write.xlsx(eTabcnFTX, paste0(saveImg, "TX_cluster_content_full_labels_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=F)

	write.xlsx(eTabcnFWT, paste0(saveImg, "WT_cluster_content_full_labels_ALL_CLUSTERS_exp_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=F)

	########################################################################################################################################################################

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
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition_ALL_CLUSTERS_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(df, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition_ALL_CLUSTERS_", expN, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	print(paste0("######## SEARCHING FOR ", cellLab[cT], " IN THE WHOLE DATASET ########"))

	# and subset the original single cell experiment using the name of the cellType cells.
	cTSubSet <- clusteredSce[, cTPos]

	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05,]

	# alternative to the filter above
	#	cTSubSet <- cTSubSet[rowSums(counts(cTSubSet) > 0) > 10,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)
	
	# exporting logcounts for Elena's violins
	subWT <- cTSubSet[, which(cTSubSet$condition=="WT")]
	subTX <- cTSubSet[, which(cTSubSet$condition=="3xTG")]

	logCsWT <- cbind.data.frame(genes=rownames(subWT), as.data.frame(logcounts(subWT)))
	logCsTX <- cbind.data.frame(genes=rownames(subTX), as.data.frame(logcounts(subTX)))

	write.xlsx(logCsWT, paste0(saveIn, "logcounts_ALL_CLUSTERS_WT_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=T)
	write.xlsx(logCsTX, paste0(saveIn, "logcounts_ALL_CLUSTERS_3xTG_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=T)

	# exporting logcounts by condition for Elena's violins
	subWT <- cTSubSet[, which(cTSubSet$condition=="WT")]
	subTX <- cTSubSet[, which(cTSubSet$condition=="3xTG")]
	
	
	# loop through all the clusters containing cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		cTSubSetAllOthers <- cTSubSet
	
		print(paste("cluster position ", cS))
	
		# subset SCE to get only cS cluster cells
		cTSubSetAllOthers$whatCluster[cTSubSetAllOthers$slm_0.5 %in% clustS[cS]] <- "thisOne"
		cTSubSetAllOthers$whatCluster[cTSubSetAllOthers$slm_0.5 %in% clustS[-cS]] <- "allOther"
		
		print(paste0("######## COMPARING ", cellLab[cT], " FROM CLUSTER ", clustS[cS], " ########"))

		cat("using findMarkers\n")

		cTSce <- list()
		conD <- unique(cTSubSetAllOthers$whatCluster)
		# then get CTRL cellType
		cTSce[[conD[1]]] <- cTSubSetAllOthers[, which(cTSubSetAllOthers$whatCluster==conD[1])]
		# and EAE
		cTSce[[conD[2]]] <- cTSubSetAllOthers[, which(cTSubSetAllOthers$whatCluster==conD[2])]
		
		if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {

			# set CTRL as a reference, to compute fold change
			cTSubSetAllOthers$whatCluster <- relevel(as.factor(cTSubSetAllOthers$whatCluster), ref="allOther")
		
			# export logcounts
			write.xlsx(cbind.data.frame(rownames(logcounts(cTSubSetAllOthers[, cTSubSetAllOthers$slm_0.5 %in% clustS[cS]])), logcounts(cTSubSetAllOthers[, cTSubSetAllOthers$slm_0.5 %in% clustS[cS]])), paste0(saveUpDown, "logcounts_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, "_cluster_", clustS[cS], ".xlsx"))
			write.xlsx(cbind.data.frame(rownames(logcounts(cTSubSetAllOthers[, cTSubSetAllOthers$slm_0.5 %in% clustS[-cS]])), logcounts(cTSubSetAllOthers[, cTSubSetAllOthers$slm_0.5 %in% clustS[-cS]])), paste0(saveUpDown, "logcounts_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, "_allOther_", paste(clustS[-cS], collapse="_"), ".xlsx"))
			
			# FINDMARKERS for DEGS
			markersG <- findMarkers(cTSubSetAllOthers, groups=cTSubSetAllOthers$whatCluster, pval.type="any", direction="any")
			
			cC <- "thisOne"
			chosen <- markersG[[cC]]

			# RELAXING THE THRESHOLD, is not done for meninges
			selfMade <- chosen[chosen$FDR < 0.05, ]
			#selfMade <- chosen[chosen$FDR < 0.1, ]

			# and export data
			finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

			# sort by logfoldchange
			sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

			# and export
#			write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_", cellLab[cT], "_CELLS_", expN, "_cluster_", clustS[cS], "_versus_all_others.csv"), col.names=T, row.names=F, sep="\t")
			write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, "_cluster_", clustS[cS], "_versus_all_others.xlsx"))
		}
		
		########################################################### WT versus 3x inside same cluster
		
		# COMPARING WT VERSUS 3XTG
		cTSubSetWTTX <- cTSubSet[, cTSubSet$slm_0.5==clustS[cS]]

		cTSce <- list()
		conD <- unique(clusteredSce$condition)
		# then get WT cellType
		cTSce[[conD[1]]] <- cTSubSetWTTX[, which(cTSubSetWTTX$condition==conD[1])]
		# and 3xTG
		cTSce[[conD[2]]] <- cTSubSetWTTX[, which(cTSubSetWTTX$condition==conD[2])]
		
		if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
		
			# set WT as a reference, to compute fold change
			cTSubSetWTTX$condition <- relevel(as.factor(cTSubSetWTTX$condition), ref="WT")

			print("finding markers")
		
			# FINDMARKERS for DEGS
			markersG <- findMarkers(cTSubSetWTTX, groups=cTSubSetWTTX$condition, pval.type="any", direction="any")

			cC <- "3xTG"
			chosen <- markersG[[cC]]

			# RELAXING THE THRESHOLD, is not done for meninges
			selfMade <- chosen[chosen$FDR < 0.05, ]
			#selfMade <- chosen[chosen$FDR < 0.1, ]

			# and export data
			finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

			# sort by logfoldchange
			sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

			# and export
	#		write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, ".csv"), col.names=T, row.names=F, sep="\t")
			write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_cluster_", clustS[cS], "_", cellLab[cT], "_CELLS_WT_versus_3xTG_", expN, ".xlsx"), overwrite=TRUE)
		}
	}
	
	########################################################### WT versus 3x all clusters

	cTSce <- list()
	conD <- unique(clusteredSce$condition)
	# then get WT cellType
	cTSce[[conD[1]]] <- cTSubSet[, which(cTSubSet$condition==conD[1])]
	# and 3xTG
	cTSce[[conD[2]]] <- cTSubSet[, which(cTSubSet$condition==conD[2])]

	if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
	
		# set WT as a reference, to compute fold change
		cTSubSet$condition <- relevel(as.factor(cTSubSet$condition), ref="WT")
		
		print("finding markers")
	
		# FINDMARKERS for DEGS
		markersG <- findMarkers(cTSubSet, groups=cTSubSet$condition, pval.type="any", direction="any")

		cC <- "3xTG"
		chosen <- markersG[[cC]]

		# RELAXING THE THRESHOLD, is not done for meninges
		selfMade <- chosen[chosen$FDR < 0.05, ]
		#selfMade <- chosen[chosen$FDR < 0.1, ]

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# sort by logfoldchange
		sortedFinal <- finalGenes[order(finalGenes$logFC, decreasing=T), ]

		# and export
#		write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_", cellLab[cT], "_CELLS_", expN, ".csv"), col.names=T, row.names=F, sep="\t")
		write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ALL_CLUSTERS_cluster_", cellLab[cT], "_CELLS_WT_vesus_3xTG_", expN, ".xlsx"), overwrite=TRUE)
	}
}
