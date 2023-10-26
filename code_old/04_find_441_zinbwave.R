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

# import plotting functions
source("00_functions.R")

# initialise label
lbl <- 1

# set path
rawPath <- "/home/gabriele/Desktop/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"

cellType <- c("Neutrophils", "CD4", "Tgd")
cellLab <-  c("Neutrophils", "CD4", "Tgd")

for (cT in seq(1, length(cellType), by=1)) {

#	cT <- 3

	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
	print("")

	# then do everything
	saveIn <- paste0(rawPath, "results_ribosomal_441/")
	saveImg <- paste0(saveIn, "figures/analysis_", cellLab[cT], "/")
	saveUpDown <- paste0(saveIn, "enrichment_", cellLab[cT], "/UpDown/")
	saveRanks <- paste0(saveIn, "enrichment_", cellLab[cT], "/ranks/")
	saveOnOff <- paste0(saveImg, "onOff/")
	saveExpression <- paste0(saveImg, "violins/")
	saveHeat <- paste0(saveImg, "heatmaps/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))
	ifelse(dir.exists(saveUpDown), TRUE, dir.create(saveUpDown, recursive=T))
	ifelse(dir.exists(saveRanks), TRUE, dir.create(saveRanks, recursive=T))
	ifelse(dir.exists(saveOnOff), TRUE, dir.create(saveOnOff, recursive=T))
	ifelse(dir.exists(saveExpression), TRUE, dir.create(saveExpression, recursive=T))
	ifelse(dir.exists(saveHeat), TRUE, dir.create(saveHeat, recursive=T))

	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")

	# NOTE: we are using QCed data file which logcounts were normalised using sizeFactors
	# but this data file contains also the information from the SLM clustering
	# performed using the zinbwaved data, in order to have a good clustering but also all
	# the genes, instead of the best 1000 used for zinbwave
	zinHV <- 1000
	clusteredSce <- readRDS(paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_20_top_", zinHV, "_withRibosomes.Rds"))
	
	# plotting up-regulated tSNEs, for specific genes
	zbWaved <- readRDS(paste0(saveIn, "slm_clustered_data_", expN, "_K_20_top_", zinHV, "_withRibosomes.Rds"))
	# first create seurat object
	genesToPlot <- c("Il2ra", "Il2rb", "Il2rg", "Il17a", "Il17f", "Il17re", "Il17ra", "Il17rc", "Il17rb", "Il17rd", "Il17d", "Il17b", "Il7r", "Bcl2a1d", "Arrb1", "Igf1r", "Rpl29", "Rps25", "S100a4", "Ccr2", "Cxcr6", "Nt5e", "Cd9", "Rorc", "Rora", "Il1r1", "Tmem176a", "Tmem176b", "Maf", "Il18r1", "Igf1r", "Ltb4r1")
	
	# adding TSNE to clusteredSce
	reducedDim(clusteredSce, "TSNE") <- reducedDim(zbWaved, "TSNE")
	
	sObj <- as.Seurat(clusteredSce, data="logcounts")
	for (gTP in genesToPlot) {
		if (gTP %in% rownames(clusteredSce)) {
			# prepare to plot with on-off genes
			sObj$gTP <- assays(clusteredSce)$logcounts[gTP, ]
			
			lbl <- lbl + 1
			# the plot
			p <- FeaturePlot(sObj, features=gTP, label=FALSE, reduction="TSNE", pt.size=1.5)
			
			# export
			png(paste0(saveOnOff, "04_0", lbl, "_", gTP, "_cluster.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
				plot(p)
			dev.off()
		} else {
			print(paste0("gene ", gTP, " not found!"))
		}
	}
	
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

	# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_clusters_matrix.png"), type="cairo", units="in", width=9, height=12, pointsize=12, res=300)
		pheatmap(log2(tab_fine_to_main+10), color=colorRampPalette(c(colorsPal[1], colorsPal[6]))(101))
	dev.off()

	# MARKER GENES DETECTION

	# three different approaches explained:
	# using ImmGen we are able to label raw, not normalised, cells to obtain their predicted cell type. Then we perform two different normalisations that are ZinbWave and sizeFactors. At this point two possibilities are open. The actual clustering is performed using the data normalised by Zinbwave while sizeFactors are not considered.
	# first option: using only cellType cells that fall into the clusters containing cellType cells
	# second option: using all cells that fall into the clusters containing cellType cells
	# third option: using all cells that were labelled by ImmGen

	# AB info is also available to validate the cells.

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
	png(paste0(saveImg, "04_0", lbl, "_", conD[2], "_cluster_composition.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	p <- plot_membership(df, cond=conD[1], col=colorsPal[6])
	# update label
	lbl <- lbl + 1
	png(paste0(saveImg, "04_0", lbl, "_", conD[1], "_cluster_composition.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()

	######################################################################################################################################################################################################################################################################################################################################################################################################################

	# declare some variables to store the markers of each cluster for each option
	optOneGenes <- list()
	optTwoGenes <- list()
	optThreeGenes <- list()

	# set threshold for FDR
	thresholD <- 0.1

	# set a personal threshold determining the minimum number of cells to compare between two conditions
	personalT <- 10

	# loop through all the clusters containing cellType
	for (cS in clustS) {

		# cS <- 10

		# OPTION ONE: findMarkers ON ALL THE IMMGEN LABELLED CELLS OF A SPECIFIC CELLTYPE IN EACH CLUSTER #

		# subset SCE to get only cS cluster cells
		optOne <- clusteredSce[, clusteredSce$slm_0.5==cS]

		print("#################################################################################")
		print(paste0("analysing cluster ", cS, " cells content: "))

		print(table(optOne$fine_to_main))

		print(paste0("######## SEARCHING FOR ", cellLab[cT], " IN CLUSTER ", cS, " ########"))

		# now find cellType only in cS
		cTPos <- grep(cellType[cT], optOne$pruned_fine)

		# and subset the original single cell experiment using the name of the cellType cells
		cTSubSet <- optOne[, cTPos]
		
		if (ncol(cTSubSet) > 0) {
			# run option one on cTSubSet
			optOneGenes <- optionOne(cTSubSet, expN, saveUpDown, personalT, conD, cellLab[cT], thresholD)
			
			# check if everything went fine
			if (length(optOneGenes) == 0) {
				print("no genes, no enrichment!")
			} else if (length(optOneGenes) > 1){
				
				# run enrichment on cTSubset
				ranksVals <- rankMyData(cTSubSet)
				# export rankings to perform later enrichment
				rTab <- cbind.data.frame(gene=names(ranksVals), deseqStat=ranksVals)
				write.table(rTab, paste0(saveRanks, "rankings_genes_for_", cellLab[cT], "_cluster_", cS, ".csv"), col.names=T, row.names=F, sep="\t")
				
				lbl <- lbl + 1
				# heatmapping
				png(paste0(saveHeat, "04_0", lbl, "_up_down_genes_heatmap_", cellLab[cT], "_cluster_", cS, ".png"), type="cairo", units="in", width=15, height=12, pointsize=12, res=300)
					scater::plotHeatmap(cTSubSet, features=optOneGenes, exprs_values="logcounts", center=TRUE, cluster_cols=FALSE, colour_columns_by=c("condition"), show_colnames=FALSE,color = rev(RColorBrewer::brewer.pal(11,name = "RdBu")))
				dev.off()
			} else {
				print("only one gene found, no heatmap!")
			}
		} else {
			print(paste0("not enough cells to compare the two conditions for ", cellLab[cT], " in cluster ", cS))
		}

		#############################################################################################

		# OPTION TWO: ALL THE CELLS IN THE CLUSTERS CONTAINING cellType CELLS

		print(paste0("######## SEARCHING FOR ALL THE CELLS IN CLUSTER ", cS, " ########"))

		# get cells in cluster cS
		cTSubSet <- clusteredSce[, clusteredSce$slm_0.5==cS]
		cTSubSet$fine_to_main[is.na(cTSubSet$fine_to_main)] <- "unknown"

		if (ncol(cTSubSet) > 0) {			
			# in case of 10x, then run this line since no AB are provided
			optTwoGenes <- optionTwo(cTSubSet, personalT, conD, cellLab[cT], thresholD)
			
			# check if everything went fine
			if (length(optTwoGenes) == 0) {
				print("no genes, no enrichment!")
			} else if (length(optTwoGenes) > 1) {
				
				lbl <- lbl + 1
				# heatmapping
				png(paste0(saveHeat, "04_0", lbl, "_up_down_genes_heatmap_ALL_cluster_", cS, ".png"), type="cairo", units="in", width=15, height=12, pointsize=12, res=300)
					scater::plotHeatmap(cTSubSet, features=optTwoGenes, exprs_values="logcounts", center=TRUE, cluster_cols=FALSE, colour_columns_by=c("condition"), show_colnames=FALSE,color = rev(RColorBrewer::brewer.pal(11,name = "RdBu")))
				dev.off()
			} else {
				print("only one gene found, no heatmap!")
			}
			
			# no GSEA here, as decided with Elena
		} else {
			print(paste0("not enough cells to compare the two conditions for cluster ", cS))
		}
	}

	#############################################################################################

	# OPTION THREE: ONLY cellType CELLS IN THE FULL DATASET

	print(paste0("######## SEARCHING FOR ", cellType[cT], " IN THE WHOLE DATASET ########"))

	# search for cellType
	cTPos <- grep(cellType[cT], clusteredSce$pruned_fine)

	# and subset the original single cell experiment using the name of the cellType cells.
	cTSubSet <- clusteredSce[, cTPos]
	
	# loop through the genes of interest
	for (gTP in genesToPlot) {
		if (gTP %in% rownames(clusteredSce)) {
			# plot the logcounts for cellType cells only, as violins, to compare the two conditions
			lbl <- lbl + 1
			p <- plotExpression(cTSubSet, x="condition", colour_by="condition", features=gTP, show_median=T) + stat_compare_means(method = "wilcox.test", label.x=1.5)
			png(paste0(saveExpression, "04_0", lbl, "comparing_", gTP, "_by_condition.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
				plot(p)
			dev.off()
		} else {
			print(paste0("gene ", gTP, " not found!"))
		}
	}

	if (ncol(cTSubSet) > 0) {
		# run option three on cTSubSet
		optThreeGenes <- optionThree(cTSubSet, expN, saveUpDown, personalT, conD, cellLab[cT], thresholD)

		################################################################### ADD SEURAT ###################################################################
		# load libraries for seurat
		library("MAST")
		library("DESeq2")
		
		# create seurat object out of the cTSubSet
		seuratSubSet <- as.Seurat(cTSubSet, data="logcounts")
		# add cell levels to seurat object
		Idents(seuratSubSet) <- cTSubSet$condition		
		
		# then perform some FindMarkers
		seuratMarkers <- list()
		seuratMarkers <- fMSeurat(seuratSubSet, seuratMarkers, thresholD, saveUpDown, expN, cellLab[cT])
		################################################################### END SEURAT ###################################################################
		
		################################################################### ADD SCDD ###################################################################
		
		library("scDD")
		cTPos <- grep(cellType[cT], clusteredSce$pruned_fine)
		# and subset the original single cell experiment using the name of the cellType cells.
		cTSubSet <- clusteredSce[, cTPos]
	
		paR <- list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)
		cTSubSet <- preprocess(cTSubSet, zero.thresh=0, scran_norm=TRUE)

		scDDegs <- scDD(cTSubSet, prior_param=paR, testZeroes=TRUE)

		################################################################### END SCDD ###################################################################
		
		# check if everything went fine
		if (length(optThreeGenes) == 0) {
			print("no genes, no enrichment!")
		} else if (length(optThreeGenes) > 1) {
			
			# now perform GSEA on cTSubSet rank file names
			lbl <- lbl + 1
			# run enrichment on cTSubset
			ranksVals <- rankMyData(cTSubSet)
			# export rankings to perform later enrichment
			rTab <- cbind.data.frame(gene=names(ranksVals), deseqStat=ranksVals)
			write.table(rTab, paste0(saveRanks, "rankings_genes_for_ALL_", cellLab[cT], "_CELLS.csv"), col.names=T, row.names=F, sep="\t")
			
			lbl <- lbl + 1
			# heatmapping
			png(paste0(saveHeat, "04_0", lbl, "_up_down_genes_heatmap_ONLY_", cellLab[cT], ".png"), type="cairo", units="in", width=15, height=12, pointsize=12, res=300)
				scater::plotHeatmap(cTSubSet, features=optThreeGenes, exprs_values="logcounts", center=TRUE, cluster_cols=FALSE, colour_columns_by=c("condition"), show_colnames=FALSE,color = rev(RColorBrewer::brewer.pal(11,name = "RdBu")))
			dev.off()
		} else {
				print("only one gene found, no heatmap!")
		}
	} else {
		print(paste0("not enough cells to compare the two conditions"))
	}
	
	############################################### PERFORMING THE RANDOM ANALYSIS
		
	genesList <- list()
	for (rnd in seq(1, 100, by=1)) {	
		
		# performing random sampling on 29 WT cells
		tXS <- sample(ncol(cTSubSet[, which(cTSubSet$condition==cnD[1])]), 29)
		# and 120 3xTG cells
		wTS <- sample(ncol(cTSubSet[, which(cTSubSet$condition==cnD[2])], 120)
		genesList[[rnd]] <- optionThreeRandomSampling(cTSubSet, expN, saveUpDown, personalT, conD, cellLab[cT], thresholD)
	}

	data <- as.data.frame(table(genesList))
	colnames(data) <- c("genes", "frequency")
	
	p <- ggplot(data, aes(x=genes, y=frequency)) +
		geom_segment( aes(x=genes, xend=genes, y=0, yend=frequency)) +
		geom_point( size=5, color="red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2) +
		coord_flip()
		
	lbl <- lbl + 1
	png(paste0(saveHeat, "04_0", lbl, "_genes_random_sampling_", cellLab[cT], ".png"), type="cairo", units="in", width=15, height=12, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	############################################### DONE PERFORMING THE RANDOM ANALYSIS

	print("                                                 ")
	print("                                                 ")
	print("_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_")
	print("                                                 ")
	print("@@@@@@@@@              END              @@@@@@@@@")
	print("                                                 ")
	print("_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_")
}
