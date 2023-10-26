#' ---
#' title: "Quality controls"
#' author: "gabrielet"
#' output: html_document
#' date: '`r format(Sys.Date(), "%d %B, %Y")`'
#' ---
#' 
#' ```{r setup, include=FALSE}
#' knitr::sharedGenes_chunk$set(
#'   tidy = TRUE,
#'   tidy.sharedGenes = list(width.cutoff = 120),
#'   message = FALSE,
#'   warning = FALSE,
#'   eval = TRUE
#' )
#' ```

library("org.Mm.eg.db")
library("clusterProfiler")
library("ggplot2")
library("RColorBrewer")
library("stringr")
library("ddpcr")
library("pathview")
library("DOSE")
library("ReactomePA")
library("enrichplot")
library("openxlsx")
library("UpSetR")

# import plotting functions
source("00_functions.R")

# set path
rawPath <- "/home/gabriele/work/cbmc/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"

# 10x cellTypes
#cTypes <-  c("CD8", "Neutrophils")

cTypes <-  c("CD8")

#for (cellType in cTypes) {

	cellType <- cTypes[1]

	# then do everything
	saveIn <- paste0(rawPath, "ANALYSIS_FINAL_441/")
	saveImg <- paste0(saveIn, "figures/analysis_", cellType, "/")
	saveGO <- paste0(saveIn, "enrichment_", cellType, "/GO_plots/")
	saveReactome <- paste0(saveIn, "enrichment_", cellType, "/Reactome_plots/")
	saveCProf <- paste0(saveIn, "enrichment_", cellType, "/cProfiler/")
	saveCProfGOBP <- paste0(saveCProf, "GO_bio_process/")
	saveCProfGOMF <- paste0(saveCProf, "GO_mol_function/")
	saveCProfRCTP <- paste0(saveCProf, "REACTOME_pathway/")
	saveCProfKEGG <- paste0(saveCProf, "KEGG_pathway/")
	saveUpDown <- paste0(saveIn, "enrichment_", cellType, "/UpDown/")
	getRanks <- paste0(saveIn, "enrichment_", cellType, "/ranks/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=T))
	ifelse(dir.exists(saveGO), TRUE, dir.create(saveGO, recursive=F))
	ifelse(dir.exists(saveReactome), TRUE, dir.create(saveReactome, recursive=F))
	ifelse(dir.exists(saveCProf), TRUE, dir.create(saveCProf, recursive=F))
	ifelse(dir.exists(saveCProfGOBP), TRUE, dir.create(saveCProfGOBP, recursive=F))
	ifelse(dir.exists(saveCProfGOMF), TRUE, dir.create(saveCProfGOMF, recursive=F))
	ifelse(dir.exists(saveCProfRCTP), TRUE, dir.create(saveCProfRCTP, recursive=F))
	ifelse(dir.exists(saveCProfKEGG), TRUE, dir.create(saveCProfKEGG, recursive=F))

	# create list of pathways to pass to enrichmydata function, to properly store the generated csv files
	# ORDER IS FUNDAMENTAL HERE
	cProfPathsList <- c(saveCProfGOBP, saveCProfGOMF, saveCProfRCTP, saveCProfKEGG)

	# set palette with 12 colours for colouring plots
	colorsPal <- brewer.pal(12, "Paired")

	# get list containing both up and down regulated genes, coming from findMarkers from one of the three optOne, optTwo, optThree functions
#	filesList <- list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))]
#	filesList <- list.files(saveUpDown)[grep("EDGER", list.files(saveUpDown))][grep("csv", list.files(saveUpDown)[grep("EDGER", list.files(saveUpDown))])]
#	filesList <- list.files(saveUpDown)[grep("FindMarkers", list.files(saveUpDown))][grep("csv", list.files(saveUpDown)[grep("FindMarkers", list.files(saveUpDown))])]

	filesListTmp <- list.files(saveUpDown)[grep("FindMarkers_wt_versus_3x_", list.files(saveUpDown))][grep("cluster", list.files(saveUpDown)[grep("FindMarkers_wt_versus_3x_", list.files(saveUpDown))])]
	filesList <- filesListTmp[grep("csv", filesListTmp)]

	# initialise label count
	lbl <- 1

	# create list of all the genes, to find those that are shared across all the clusters
	sharedGenes <- list()

	# loop through the gene lists
	for (fName in filesList) {
	
#		fName <- filesList[1]

		print(paste0("analysing ", fName))
		#load data
		originalData <- read.csv(paste0(saveUpDown, fName), header=T, sep="\t", stringsAsFactors=T)
		# get identifier for this experiment
#		ID <- strsplit(fName, "_up")[[1]][1]
		ID <- strsplit(fName, paste0("_", cellType))[[1]][1]

		# check if originalData contains actual information
		if (nrow(originalData) != 0) {
			# set thresholds for enrichment analyses
			# first value is pvalue cutoff, second value is qvalue cutoff
			threshS <- c(0.2, 0.2)
			enrichMyData(originalData, cProfPathsList, ID, threshS, expN)
		} else {
			print("no genes to enrich")
		}
		
		# LOLLIPLOTTING
		
#		GOs <- read.csv(paste0(cProfPathsList[1], "cProf_", ID, "_GO_bio_process_", expN, ".csv"), header=T, sep="\t")
#		KEGGs <- read.csv(paste0(cProfPathsList[4], "FindMarkers_KEGG_441_LA_VERSIONE_DI_ELENA.csv"), header=T, sep="\t")
#		
#		# threshold for gene count
#		thresh <- c(4, 8)
#		# and limits for x axis. the lower the threshold, the longer the axis (i.e. more occurences!)
#		limS <- c(45, 30)
#		
#		for (thRS in c(1, 2)) {
#		
#			# get interesting GOs and KEGGs
#			tmpGOs <- GOs[GOs$p.adjust <= 0.05, ]
#			intGOs <- tmpGOs[tmpGOs$Count >= thresh[thRS], ]
#			
#			tmpKEGGs <- KEGGs[KEGGs$p.adjust <= 0.05, ]
#			intKEGGs <- tmpKEGGs[tmpKEGGs$Count >= thresh[thRS], ]
#			
#			# group by categories		
#			catOneGO <- length(grep("dhesion|igration|hemotaxis|ctin|hemokine", intGOs$Description))
#			catTwoGO <- length(grep("ntigen|ntigenic|mmune response", intGOs$Description))
#			catThreeGO <- length(grep("ymphocyte differentiation|T cell differentiation", intGOs$Description))
#			catFourGO <- length(grep("ytokine|nterleukin|umor necrosis factor", intGOs$Description))
#			catFiveGO <- length(grep("nflammation|nflammatory", intGOs$Description))

#			catOneKEGG <- length(grep("Infection", intKEGGs$elenaLabels))
#			catTwoKEGG <- length(grep("Neurodegenerative disease", intKEGGs$elenaLabels))
#			catThreeKEGG <- length(grep("Immune response", intKEGGs$elenaLabels))
#			catFourKEGG <- length(grep("Th17|T helper 17|IL-17", intKEGGs$Description))
#			
#			allCatSGO <- cbind.data.frame(categories=c("catOneGO", "catTwoGO", "catThreeGO", "catFourGO", "catFiveGO"), sumTot=c(catOneGO, catTwoGO, catThreeGO, catFourGO, catFiveGO), coloR =c("GO", "GO", "GO", "GO", "GO"))
#			catSortGO <- allCatSGO[order(allCatSGO$sumTot), ]
#			
##			catSortGO$categories <- with(catSortGO, reorder(categories, sumTot))
#			
#			allCatSKEGG <- cbind.data.frame(categories=c("catOneKEGG", "catTwoKEGG", "catThreeKEGG", "catFourKEGG"), sumTot=c(catOneKEGG, catTwoKEGG, catThreeKEGG, catFourKEGG), coloR=c("KEGG", "KEGG", "KEGG", "KEGG"))
#			catSortKEGG <- allCatSKEGG[order(allCatSKEGG$sumTot), ]
#			
##			catSortKEGG$categories <- with(catSortKEGG, reorder(categories, sumTot))
#			
#			catSort <- rbind.data.frame(catSortKEGG, catSortGO)
#			catSort$categories <- factor(catSort$categories, levels=unique(catSort$categories))
#			
#			# plot			
#			p <- ggplot(catSort, aes(x=categories, y=sumTot, fill=coloR)) + geom_bar(stat = "identity", width=0.45, color="black") + scale_y_continuous(limits=c(0, limS[thRS]), breaks=seq(0, limS[thRS], by=5)) + xlab("") + ylab("") + coord_flip() + theme(aspect.ratio = 1)

#				
#			png(paste0(saveImg, "05_01_BARPLOT_", expN, "_threshold_", thresh[thRS], "_RATIO_1.png"), type="cairo", units="in", width=8, height=8, pointsize=12, res=300)
#				plot(p)
#			dev.off()
#			
#			p <- ggplot(catSort, aes(x=categories, y=sumTot, fill=coloR)) + geom_bar(stat = "identity", width=0.8, color="black") + scale_y_continuous(limits=c(0, limS[thRS]), breaks=seq(0, limS[thRS], by=5)) + xlab("") + ylab("") + coord_flip() + theme(aspect.ratio = 1)

#				
#			png(paste0(saveImg, "05_01_BARPLOT_", expN, "_threshold_", thresh[thRS], "_RATIO_1_BARS_BIG.png"), type="cairo", units="in", width=8, height=8, pointsize=12, res=300)
#				plot(p)
#			dev.off()
#					
#			p <- ggplot(catSort, aes(x=categories, y=sumTot, fill=coloR)) + geom_bar(stat = "identity", width=0.45, color="black") + scale_y_continuous(limits=c(0, limS[thRS]), breaks=seq(0, limS[thRS], by=5)) + xlab("") + ylab("") + coord_flip() + theme(aspect.ratio = 0.2)

#				
#			png(paste0(saveImg, "05_01_BARPLOT_", expN, "_threshold_", thresh[thRS], "_RATIO_02.png"), type="cairo", units="in", width=8, height=8, pointsize=12, res=300)
#				plot(p)
#			dev.off()
#			
#		}
#############################################		# get cluster name for upset plotting
#############################################		fN <- strsplit(fName, "_exp")[[1]][1]
#############################################		
#############################################		# and store the genes, if there are some
#############################################		if (length(originalData$gName) > 0) {
#############################################			sharedGenes[[fN]] <- originalData$gName
#############################################		}
	}

#############################################	# upset plot with the intersection of the up_and_down genes

#############################################	# are there shared genes?
#############################################	intersection <- unique(unlist(sharedGenes))

#############################################	# export the full list with shared genes
#############################################	write.table(intersection, paste0(saveImg, "shared_up_and_down_genes.csv"), col.names=F, row.names=F, sep="\t")

#############################################	# to perform a comparison al least two lists are required
#############################################	if (length(names(sharedGenes)) > 1) {

#############################################		# loop through the different sets of genes
#############################################		for (lsT in names(sharedGenes)) {

#############################################			# if at least an intersection is found
#############################################			if (length(intersection) != 0) {
#############################################			
#############################################				uP <- upset(fromList(sharedGenes), order.by="freq", nsets=length(names(sharedGenes)), decreasing=T, sets=names(sharedGenes), keep.order=T, text.scale=1.5, sets.x.label="Set size", mainbar.y.label="Shared proteins")	
#############################################				
#############################################				# update label
#############################################				lbl <- lbl + 1
#############################################				png(paste0(saveImg, "04_0", lbl, "_upset_", lsT), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
#############################################					print(uP)
#############################################				dev.off()

#############################################				# plotting heatmap
#############################################				geneS <- vector()
#############################################				clustS <- vector()
#############################################				belonG <- vector()
#############################################				# loop through the cols
#############################################				for (nC in seq(1, ncol(uP$New_data), by=1)) {
#############################################					# loop through the rows
#############################################					for (nR in seq(1, nrow(uP$New_data), by=1)) {
#############################################						# get cluster name
#############################################						clustS <- append(clustS, names(uP$New_data)[nC])
#############################################						# get gene name
#############################################						geneS <- append(geneS, rownames(uP$New_data)[nR])
#############################################						# get belonging 
#############################################						bG <- uP$New_data[nR, nC]
#############################################						if (bG == 0) {
#############################################							belonG <- append(belonG, "not found")
#############################################						} else {
#############################################							belonG <- append(belonG, "found")
#############################################						}
#############################################					}
#############################################				}
#############################################				# put everything together
#############################################				dataG <- cbind.data.frame(geneS, clustS, belonG)

#############################################				# and print 
#############################################				p <- ggplot(dataG, aes(clustS, geneS, fill=factor(belonG))) + 
#############################################					geom_tile(color="black") +
#############################################					xlab("cluster") +
#############################################					ylab("gene name") +
#############################################					guides(fill=guide_legend(title="Gene in cluster is: "), position="bottom") +
#############################################					scale_fill_manual(values=c("not found"=colorsPal[1], "found"=colorsPal[6]))

#############################################				# update label
#############################################				lbl <- lbl + 1
#############################################				png(paste0(saveImg, "04_0", lbl, "_heatmap_", lsT), type="cairo", units="in", width=6, height=32, pointsize=12, res=300)
#############################################					plot(p)
#############################################				dev.off()
#############################################			}
#############################################		}
#############################################	} else {
#############################################		print("the number of list is not sufficient to perform comparisons!")
#############################################	}

#############################################	# get list containing both up and down regulated genes
#############################################	filesList <- list.files(getRanks)

#############################################	# set threshold for DESeq2 statistics
#############################################	deSeqStatThresh <- 2

#############################################	# loop through the gene lists
#############################################	for (fName in filesList) {

#############################################		print(paste0("analysing ", fName))
#############################################		
#############################################		lbl <- 1
#############################################		
#############################################		# load data
#############################################		originalData <- read.csv(paste0(getRanks, fName), header=T, sep="\t", stringsAsFactors=T)
#############################################		
#############################################		# filter originalData using fold change
#############################################		originalData <- originalData[which(abs(originalData$deseqStat) >= deSeqStatThresh), ]
#############################################		
#############################################		# check if we are working with a cluster
#############################################		if (length(grep("cluster", fName)) != 0) {
#############################################			# get cluster name
#############################################			cS <- strsplit(strsplit(fName, "cluster_")[[1]][2], "\\.")[[1]][1]
#############################################			# check if originalData contains actual information
#############################################			if (nrow(originalData) != 0) {
#############################################				lbl <- goPlotting(originalData, saveGO, saveCProf, lbl, fName, cS)
#############################################				lbl <- reactomePlotting(originalData, saveReactome, saveCProf, lbl, fName, cS)
#############################################			} else {
#############################################				print("no genes to enrich")
#############################################			}
#############################################		# or with all the cells
#############################################		} else if (length(grep("cluster", fName)) == 0) {
#############################################			# set cluster name
#############################################			cS <- "ALL_CELLS"
#############################################			# check if originalData contains actual information
#############################################			if (nrow(originalData) != 0) {
#############################################				lbl <- goPlotting(originalData, saveGO, saveCProf, lbl, fName, cS)
#############################################				lbl <- reactomePlotting(originalData, saveReactome, saveCProf, lbl, fName, cS)
#############################################			} else {
#############################################				print("no genes to enrich")
#############################################			}
#############################################		}
#############################################	}

	# shut all the residual open graphics down
	graphics.off()
#}
