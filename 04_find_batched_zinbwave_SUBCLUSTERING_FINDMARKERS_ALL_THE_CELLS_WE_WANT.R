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
 
# select cell types
cellType <- c("T.8|T.CD8", "Neutrophils", "Tgd")#, "Monocytes", "Macrophages", "Microglia", "B cells", "ILC", "NK cells", "NKT", "Mast")
cellLab <-  c("CD8", "Neutrophils", "Tgd")#, "Monocytes", "Macrophages", "Microglia", "B_cells", "ILC", "NK_cells", "NKT", "Mast")

# loop through cell type
for (cT in seq(1, length(cellLab), by=1)){

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
	print("")

	# then do everything
	saveIn <- paste0(rawPath, "ANALYSIS_FINAL_441/")
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
	
	cTSubSet <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, "_", cellLab[cT], "_ONLY.Rds"))

	# find the clusters
	fR <- table(cTSubSet$slm_0.3)
	clustS <- names(fR[which(fR!=0)])

	# get label ready
	cTSubSet$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(cTSubSet$pruned_fine, "\\("), `[[`, 1)))
	
	##################### MAIN labels	
	
	# plotting heatmap with cell type
	tab_fine_to_main <- table(Assigned=cTSubSet$fine_to_main, Cluster=cTSubSet$slm_0.3, Condition=cTSubSet$condition)

	# get 3x table ready to be exported
	eTabTX <- as.data.frame(tab_fine_to_main[, , "3xTG"])
	eTabRNTX <- cbind.data.frame(cellType=rownames(eTabTX), eTabTX)
	eTabCNTX <- rbind.data.frame(c("Cluster", "Freq"), eTabRNTX)

	# same for WT
	eTabWT <- as.data.frame(tab_fine_to_main[, , "WT"])
	eTabRNWT <- cbind.data.frame(cellType=rownames(eTabWT), eTabWT)
	eTabCNWT <- rbind.data.frame(c("Cluster", "Freq"), eTabRNWT)

	# export both tables
	write.xlsx(eTabCNTX, paste0(saveImg, "TX_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], ".xlsx"), colNames=F, overwrite=T)

	write.xlsx(eTabCNWT, paste0(saveImg, "WT_cluster_SUBCLUSTERED_content_exp_", expN, "_", cellLab[cT], ".xlsx"), colNames=F, overwrite=T)

	##################### FULL labels

	tab_full_labels <- table(Assigned=cTSubSet$pruned_fine, Cluster=cTSubSet$slm_0.3, Condition=cTSubSet$condition)

	# load objs
	eTabcnFTX <- NULL
	eTabcnFWT <- NULL
	if (length(unique(cTSubSet$pruned_fine)) > 1 && length(unique(cTSubSet$slm_0.3)) > 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_full_labels[, , "3xTG"])
		eTabrnFTX <- cbind.data.frame(cellType=rownames(etabFTX), etabFTX)
		eTabcnFTX <- eTabrnFTX[, c(2,3,4)]
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_full_labels[, , "WT"])
		eTabrnFWT <- cbind.data.frame(cellType=rownames(etabFWT), etabFWT)
		eTabcnFWT <- eTabrnFWT[, c(2,3,4)]
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$pruned_fine)) == 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_full_labels[, , "3xTG"])
		eTabcnFTX <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFTX)), rownames(etabFTX), etabFTX)
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_full_labels[, , "WT"])
		eTabcnFWT <- cbind.data.frame(rep(unique(cTSubSet$pruned_fine), length(etabFWT)), rownames(etabFWT), etabFWT)
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	} else if (length(unique(cTSubSet$slm_0.3)) == 1) {
		# get 3xTG full table ready to be exported
		etabFTX <- as.data.frame(tab_full_labels[, , "3xTG"])
		eTabcnFTX <- cbind.data.frame(rownames(etabFTX), rep(unique(cTSubSet$slm_0.3), length(etabFTX)), etabFTX)
		colnames(eTabcnFTX) <- c("CellType", "Cluster", "Freq")

		# and WT
		etabFWT <- as.data.frame(tab_full_labels[, , "WT"])
		eTabcnFWT <- cbind.data.frame(rownames(etabFWT), rep(unique(cTSubSet$slm_0.3), length(etabFWT)), etabFWT)
		colnames(eTabcnFWT) <- c("CellType", "Cluster", "Freq")
	}
	
	
	# export the full table
	write.xlsx(eTabcnFTX, paste0(saveImg, "TX_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], ".xlsx"), colNames=T, overwrite=T)

	write.xlsx(eTabcnFWT, paste0(saveImg, "WT_cluster_content_SUBCLUSTERED_full_labels_exp_", expN, "_", cellLab[cT], ".xlsx"), colNames=T, overwrite=T)

	####################### CLEAN AND NORMALISE #######################
	
	# remove genes with too many zeros
	cTSubSet <- cTSubSet[rowMeans(assays(cTSubSet)$counts > 0) > 0.05,]

	# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
	# compute sizeFactors
	cTSubSet <- computeSumFactors(cTSubSet, min.mean=0.5)
	cTSubSet <- logNormCounts(cTSubSet)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################

	superfull_table <- assay(cTSubSet, "logcounts")
	colnames(superfull_table) <- paste(colnames(cTSubSet), cTSubSet$condition, cTSubSet$region, cTSubSet$pruned_fine, "cluster", cTSubSet$slm_0.3, sep="_")
	with_row_names <- cbind.data.frame(rownames(superfull_table), superfull_table)
	colnames(with_row_names) <- c("genes", colnames(superfull_table))
	
	write.xlsx(with_row_names, paste0(saveIn, "logcounts_all_the_selected_cells_version_barbara_elena_eleonora_only_SUBCLUSTERED_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T)
	
	################################################# EXPORT THE NEW TABLE ELENA ASKED #################################################
	
	# exporting logcounts for Elena's violins BY CONDITION
	subWT <- cTSubSet[, which(cTSubSet$condition=="WT")]
	subTX <- cTSubSet[, which(cTSubSet$condition=="3xTG")]

	logCsWT <- cbind.data.frame(genes=rownames(subWT), as.data.frame(logcounts(subWT)))
	logCsTX <- cbind.data.frame(genes=rownames(subTX), as.data.frame(logcounts(subTX)))

	write.xlsx(logCsWT, paste0(saveIn, "logcounts_SUBCLUSTERED_WT_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=T)
	write.xlsx(logCsTX, paste0(saveIn, "logcounts_SUBCLUSTERED_3xTG_", expN, "_", cellLab[cT], ".xlsx"), overwrite=T, colNames=T)
	
	####################### TEST 1 IN EACH CLUSTER COMPARE STUFF #######################
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		# PERFORM WT VS 3X
	
		print(paste("cluster position ", clustS[cS]))
	
		# subset SCE to get only cS cluster cells
		aClust <- cTSubSet[, cTSubSet$slm_0.3==clustS[cS]]
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ", clustS[cS], " ########"))
		cTSce <- list()
		conD <- unique(aClust$condition)
		# then get CTRL cellType
		cTSce[[conD[1]]] <- aClust[, which(aClust$condition==conD[1])]
		# and EAE
		cTSce[[conD[2]]] <- aClust[, which(aClust$condition==conD[2])]
		
		if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {

			# set reference condition
			aClust$condition <- relevel(as.factor(aClust$condition), ref="WT")

			# export logcounts
#			write.xlsx(cbind.data.frame(rownames(aClust), logcounts(aClust)), paste0(saveUpDown, "logcounts", cellType, "_CELLS_", expN, "_cluster_", clustS[cS], "_", cellLab[cT], ".xlsx"), overwrite=T)
			
			cat("finding markers\n")
			
			# FINDMARKERS for DEGS
			markersG <- findMarkers(aClust, groups=aClust$condition, pval.type="any", direction="any")
	
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
#			write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_", expN, "_cluster_", clustS[cS], ".csv"), row.names=F, quote=F, sep="\t")
			write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_SUBCLUSTERED_", cellLab[cT], "_CELLS_", expN, "_cluster_wt_versus_3xtg", clustS[cS], ".xlsx"), overwrite=T)
		}
	}
	
	####################### TEST 2 COMPARE A CLUSTER AGAINST ALL OTHERS #######################
	
	print("####################### TEST 2 #######################")
	
	# get cTSubSet, i.e. all the dataset
	cTSubSetTestTwo <- cTSubSet
	
	# loop through all the clusters containing a specific cellType
	for (cS in seq(1, length(clustS), by=1)) {
	
		print(paste("cluster position ", clustS[cS]))
	
		# subset SCE to get only cS cluster cells
		# this is for comparing one cluster versus all the other clusters
		cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.3 %in% clustS[cS]] <- "thisOne"
		cTSubSetTestTwo$whatCluster[cTSubSet$slm_0.3 %in% clustS[-cS]] <- "allOther"

		# otherwise we may compare, inside a cluster, WT versus AD
		
		print(paste0("######## COMPARING ", cellType[cT], " INSIDE CLUSTER ", clustS[cS], " ########"))

		cTSce <- list()
		conD <- unique(cTSubSetTestTwo$whatCluster)
		
		if (length(conD)>1) {
		
			# then get CTRL cellType
			cTSce[[conD[1]]] <- cTSubSetTestTwo[, which(cTSubSetTestTwo$whatCluster==conD[1])]
			# and EAE
			cTSce[[conD[2]]] <- cTSubSetTestTwo[, which(cTSubSetTestTwo$whatCluster==conD[2])]
			
			if (ncol(cTSce[[conD[1]]]) >= 10 & ncol(cTSce[[conD[2]]]) >= 10) {
			
				# set ALLOTHER as a reference, to compute fold change
				cTSubSetTestTwo$whatCluster <- relevel(as.factor(cTSubSetTestTwo$whatCluster), ref="allOther")
						
				cat("finding markers\n")
				
				# FINDMARKERS for DEGS
				markersG <- findMarkers(cTSubSetTestTwo, groups=cTSubSetTestTwo$whatCluster, pval.type="any", direction="any")
		
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
#					write.table(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_SUBLCUSTERED_", cellLab[cT], "_CELLS_", expN, "_cluster_", clustS[cS], ".csv"), row.names=F, quote=F, sep="\t")
				write.xlsx(sortedFinal, paste0(saveUpDown, "up_down_FindMarkers_ONE_VERSUS_ALL_OTHER_SUBLCUSTERED_", cellLab[cT], "_CELLS_", expN, "_cluster_", clustS[cS], ".xlsx"), overwrite=T)
			}
		} else {
			print("not enough conditions to compare")
		}
	}
}
