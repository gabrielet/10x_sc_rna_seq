#' ---
#' title: "Quality controls"
#' author: "gabrielet"
#' output: html_document
#' date: '`r format(Sys.Date(), "%d %B, %Y")`'
#' ---
#' 
#' ```{r setup, include=FALSE}
#' knitr::opts_chunk$set(
#'   tidy = TRUE,
#'   tidy.opts = list(width.cutoff = 120),
#'   message = FALSE,
#'   warning = FALSE,
#'   eval = TRUE
#' )
#' ```

library("igraph")
library("ggplot2")
library("RColorBrewer")
library("STRINGdb")
library("ggrepel")
library("org.Mm.eg.db")
library("clusterProfiler")
library("ddpcr")
library("openxlsx")

firstUp <- function(strVct) {
	
	reS <- vector()
	for (stR in strVct) {
		stR <- tolower(stR)
		substr(stR, 1, 1) <- toupper(substr(stR, 1, 1))
		reS <- append(reS, stR)
	}
	return(reS)
}

# import plotting functions
source("00_functions.R")

# cellTypes
cellTypes <- c("CD4", "Tgd", "Neutrophils")

# set path
rawPath <- "/home/gabriele/Desktop/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"

# loop through cellTypes
#for (cType in cellTypes){

	# analyse gamma delta only
	cType <- cellTypes[2]

	# then do everything
	saveIn <- paste0(rawPath, "ANALYSIS_FINAL_441/")
	saveNets <- paste0(saveIn, "enrichment_", cType, "/networks/")
	saveUpDown <- paste0(saveIn, "enrichment_", cType, "/UpDown/")

	# check if saveIn directory exist, otherwise create it
	ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
	ifelse(dir.exists(saveNets), TRUE, dir.create(saveNets, recursive=F))

	# set palette with 12 colours for colouring plots
	coloursPal <- brewer.pal(12, "Paired")

	# get list containing both up and down regulated genes but only on the full dataset, not by single cluster
#	fName <- list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))][grep("ONLY", list.files(saveUpDown)[grep("up_and_down.csv", list.files(saveUpDown))])]

	# up_down_441_EDGER_Tgd_CELLS.csv
#	fName <- list.files(saveUpDown)[grep(paste0("CELLS_", expN, "_RELAXED_THRESHOLD.csv"), list.files(saveUpDown))]
	fName <- list.files(saveUpDown)[grep(paste0("CELLS_", expN, ".csv"), list.files(saveUpDown))]

	# load STRING to SYMBOL mapping file
	# this file is obtained by the info file downloaded from STRING. from this file extract the first two columns to avoid importing issues
	mappingNames <- read.table("/home/gabriele/cbmc/scrnaseq/scrnaseq/10090.protein.info_cols_reduced.v11.0.txt", header=T)

	#load data
	originalData <- read.csv(paste0(saveUpDown, fName), header=T, sep="\t", stringsAsFactors=T)
	
	# 10090 refers to Mus musculus
	stringDB <- STRINGdb$new(version="11", species=10090, score_threshold=400, input_directory=getwd())

	# mapping the genes into the STRING database to build the networks
	mappedGenes <- stringDB$map(originalData, "gName", removeUnmappedRows=F)
	
	# if mappedgenes are uppercase, make them lower with only first capital
	mappedGenes$gName <- firstUp(mappedGenes$gName)

	# found and not found genes
	notFound <- mappedGenes[is.na(mappedGenes$STRING_id), ]
	foundGenes <- mappedGenes[!is.na(mappedGenes$STRING_id), ]
	# export not found genes!
	write.table(notFound$gName, paste0(saveNets, "notFoundSTRINGGenes.csv"), row.names=FALSE, col.names=F, quote=FALSE)

	# build the network, using the geneNames (using goGroup or keggGroup geneNames should be equivalent!)
	originalNet <- igraph::simplify(stringDB$get_subnetwork(foundGenes$STRING_id), remove.multiple=TRUE, remove.loops=TRUE)
	V(originalNet)$symbol <- mappingNames$preferred_name[mappingNames$protein_external_id %in% V(originalNet)$name]
		
	# find correspondence between the corrected string names and the original names
	# get STRING ID and name
	stringCorrected <- cbind.data.frame(STRING_id=V(originalNet)$name, symbols=V(originalNet)$symbol)
	# get original names and merge with those from STRING
	combinedNames <- merge(stringCorrected, foundGenes, by="STRING_id")[, c(1, 2, 3)]
	colnames(combinedNames) <- c("stringID", "stringSymbol", "originalName")
	
	# now use these info to rename the nodes in the network, using the original names
	for (strID in seq(1, length(V(originalNet)$name), by=1)) {
		V(originalNet)$symbol[strID] <- combinedNames$originalName[which(combinedNames$stringID==V(originalNet)$name[strID])]
	}
	
	# AT THIS POINT THERE ARE SEVERAL OPTIONS TO BUILD DIFFERENT NETWORKS	
		
########	################################## analyse the connected component using logFC as node's weight
########	
########	# get the biggest connected component from the originalNet
########	connComp <- decompose(originalNet)[[1]]
########	
########	# find min logFC and compute its absolute
########	minFC <- abs(min(originalData$logFC))
########	originalNorm <- originalData
########	# transpose negative logFCs to positive, and the positive to even-more-positive
########	originalNorm$logFC <- originalData$logFC + minFC + 1
########		
########	# initialise weighted degree storage	
########	degrees <- vector()
########	# computing degree for each node in the connComp, using logFC
########	for (vrtx in V(connComp)) {
########		# get neighbours of a vertex
########		neighList <- neighbors(connComp, vrtx, mode="all")
########		# initialise degree of a vertex
########		sumUp <- 0
########		# for each neighbour of the selected vertex
########		for (ngh in neighList$symbol) {
########			# get logFC and sum this value up
########			sumUp <- sumUp + originalNorm$logFC[which(originalNorm$gName==ngh)]
########		}
########		# another vertex degree is computed
########		degrees <- append(degrees, sumUp)
########	}
########	# finally, build weighted degree table
########	degreeTable <- cbind.data.frame(geneName=V(connComp)$symbol, weightedDeg=degrees)

########	# get only nodes with degree higher to perC percentile
########	perC <- 0.95
########	highNodes <- degreeTable[which(degreeTable$weightedDeg > quantile(degreeTable$weightedDeg, perC)), ]
########	
########	# set size and symbols
########	# for those nodes NOT IN highNodes, set empty symbol
########	V(connComp)$symbol[!(V(connComp)$symbol %in% highNodes$geneName)] <- ""
########	# same for size. those IN highNodes: size = 2, others size = 0.5
########	V(connComp)$size <- 0.5
########	V(connComp)$size[V(connComp)$symbol %in% highNodes$geneName] <- 2
########	
########	# compute layout	
########	set.seed(131)
########	aR <- vcount(connComp)^2
########	layO <- layout_with_lgl(connComp, maxiter=10000, maxdelta=vcount(connComp),  area=aR, coolexp=1.5, repulserad=aR*vcount(connComp), cellsize=sqrt(sqrt(aR)), root=NULL)
########	
########	# and print network
########	png(paste0(saveNets, "04_01_degree_weighted.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
########		plot(igraph::simplify(connComp), vertex.size=V(connComp)$size, vertex.label=V(connComp)$symbol, vertex.label.cex=0.3, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5], layout=layO)
########	dev.off()
########	
########	################################## compute degree and betweenness on the connected component
########	
########	# get the biggest connected component from the originalNet
########	connComp <- decompose(originalNet)[[1]]

########	# computing Degree
########	deg <- degree(connComp, v=V(connComp), mode="all", normalized=F)
########	
########	# computing Betweenness on both networks
########	betw <- betweenness(connComp, v=V(connComp), directed=F, weights=NULL, normalized=F)

########	# obtain deg and betw names
########	names(deg) <- mappingNames$preferred_name[mappingNames$protein_external_id %in% names(deg)]
########	names(betw) <- mappingNames$preferred_name[mappingNames$protein_external_id %in% names(betw)]
########	
########	highNodes <- list()
########	# get high degree nodes, using the third quantile as threshold
########	highNodes$degree <- names(deg[which(deg > quantile(deg, perC))])
########	# get high betweenness
########	highNodes$betweenness <- names(betw[which(betw > quantile(betw, perC))])
########	
########	# find the intesection between the two 
########	interS <- intersect(highNodes$degree, highNodes$betweenness)
########	
########	# at this point, if there are some intersecting nodes
########	if (length(interS) > 0) {
########		# build data frame with the necessary info
########		toScatter <- data.frame(Degree=deg[interS], Betwenness=betw[interS])
########		
########		p <- ggplot(toScatter, aes(x=Degree, y=Betwenness)) +
########			geom_point(alpha=0.8, shape=21, color="black", fill=coloursPal[10], size=5) +
########			geom_text_repel(label=rownames(toScatter), fontface="italic", size=4, colour="black", max.overlaps=50)
########		png(paste0(saveNets , "04_01_scattering_betw_deg.png"), type="cairo", units="in", width=16, height=12, pointsize=12, res=300)
########			plot(p)
########		dev.off()
########		
########		# set size and symbols
########		# for those nodes NOT IN interS, set empty symbol
########		V(connComp)$symbol[!(V(connComp)$symbol %in% interS)] <- ""
########		# same for size. those IN interS: size = 2, others size = 0.5
########		V(connComp)$size <- 0.5
########		V(connComp)$size[V(connComp)$symbol %in% interS] <- 2
########		
########		# compute layout
########		set.seed(131)
########		aR <- vcount(connComp)^2
########		layO <- layout_with_lgl(connComp, maxiter=10000, maxdelta=vcount(connComp),  area=aR, coolexp=1.5, repulserad=aR*vcount(connComp), cellsize=sqrt(sqrt(aR)), root=NULL)
########		
########		# and print network
########		png(paste0(saveNets, "04_01_degree_original.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
########			plot(igraph::simplify(connComp), vertex.size=V(connComp)$size, vertex.label=V(connComp)$symbol, vertex.label.cex=0.3, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5], layout=layO)
########		dev.off()
########	} else {
########		print("no intersection found!")
########	}
########	
########	################################## print interleukins only, from the connected component
########	
########	# get the biggest connected component from the originalNet
########	connComp <- decompose(originalNet)[[1]]
########	
########	interLPos <- grep("Il", V(connComp)$symbol)
########	interL <- V(connComp)$symbol[interLPos]
########	
########	# set size and symbols
########	# for those nodes NOT IN interL, set empty symbol
########	V(connComp)$symbol[!(V(connComp)$symbol %in% interL)] <- ""
########	# same for size. those IN interL: size = 2, others size = 0.5
########	V(connComp)$size <- 0.5
########	V(connComp)$size[V(connComp)$symbol %in% interL] <- 2
########	
########	# compute layout
########	set.seed(131)
########	aR <- vcount(connComp)^2
########	layO <- layout_with_lgl(connComp, maxiter=10000, maxdelta=vcount(connComp),  area=aR, coolexp=1.5, repulserad=aR*vcount(connComp), cellsize=sqrt(sqrt(aR)), root=NULL)
########	
########	# and print network
########	png(paste0(saveNets, "04_01_interleukins.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
########		plot(igraph::simplify(connComp), vertex.size=V(connComp)$size, vertex.label=V(connComp)$symbol, vertex.label.cex=0.3, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5], layout=layO)
########	dev.off()
########	
########	################################## print and analyse interleukins network, i.e. with neighbours, from the connected component
########	
########	# get the biggest connected component from the originalNet
########	connComp <- decompose(originalNet)[[1]]
########	
########	# get interleukins neighbours
########	neighList <- list()
########	for (nd in seq(1, length(V(connComp)$name[interLPos]), by=1)) {
########		neighList[[nd]] <- neighbors(connComp, V(connComp)$name[interLPos][nd], mode="all")
########	}
########	
########	# build the network with all the neighbours
########	interLNet <- igraph::simplify(stringDB$get_subnetwork(names(unlist(neighList))), remove.multiple=TRUE, remove.loops=TRUE)
########	# get original names using the combinedNames table (see above)
########	for (nd in seq(1, length(V(interLNet)$name), by=1)) {
########		V(interLNet)$symbol[nd] <- combinedNames$originalName[which(combinedNames$stringID==V(interLNet)$name[nd])]
########	}
########	
########	# compute network degree to test which interleukin is more central	
########	degI <- degree(interLNet, v=V(interLNet), mode="all", normalized=F)
########	
########	# set size and symbols
########	# for those nodes NOT IN interL, set empty symbol
########	V(interLNet)$symbol[!(V(interLNet)$symbol %in% interL)] <- ""
########	# same for size. those NOT IN interL set size = 0.5
########	V(interLNet)$size <- (degI/10) + 0.5
########	V(interLNet)$size[!(V(interLNet)$symbol %in% interL)] <- 0.5
########	
########	# compute layout
########	set.seed(131)
########	aR <- vcount(interLNet)^2
########	layO <- layout_with_lgl(interLNet, maxiter=10000, maxdelta=vcount(interLNet),  area=aR, coolexp=1.5, repulserad=aR*vcount(interLNet), cellsize=sqrt(sqrt(aR)), root=NULL)
########	
########	# and print network
########	png(paste0(saveNets, "04_01_interleukins_neighs.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
########		plot(igraph::simplify(interLNet), vertex.size=V(interLNet)$size, vertex.label=V(interLNet)$symbol, vertex.label.cex=0.3, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5], layout=layO)
########	dev.off()

########	################################## find Elena's selected genes, from the connected component
########	
########	# get the biggest connected component from the originalNet
########	connComp <- decompose(originalNet)[[1]]
########	
########	# set genes list
########	geneS <- c("Il2ra", "Arrb1", "Igf1r", "Rpl29", "Rps25", "Ccr2", "Cxcr6", "Il1r1", "Maf", "Il18r1", "Il17a", "Sod1", "Cd44", "Ctse", "Ctsd", "Scart1", "Blk", "Sox13", "Nr4a1", "S100a8", "Il17rc", "Il23r", "Cd69", "Cd96")
########	
########	# get genes positions
########	gnPos <- which(V(connComp)$symbol %in% geneS == T)

########	# get genes neighbours
########	neighList <- list()
########	for (nd in seq(1, length(V(connComp)$name[gnPos]), by=1)) {
########		neighList[[nd]] <- neighbors(connComp, V(connComp)$name[gnPos][nd], mode="all")
########	}
########	
########	# build the network
########	genesNet <- igraph::simplify(stringDB$get_subnetwork(names(unlist(neighList))), remove.multiple=TRUE, remove.loops=TRUE)
########	
########	# get biggest, connected component
########	connCompRed <- decompose(genesNet)[[1]]
########	
########	# get original names
########	for (nd in seq(1, length(V(connCompRed)$name), by=1)) {
########		V(connCompRed)$symbol[nd] <- combinedNames$originalName[which(combinedNames$stringID==V(genesNet)$name[nd])]
########	}

########	# compute network degree to test which interleukin is more central	
########	degGn <- degree(connCompRed, v=V(connCompRed), mode="all", normalized=F)
########	
########	# set size and symbols
########	# for those nodes NOT IN interL, set empty symbol
########	V(connCompRed)$symbol[!(V(connCompRed)$symbol %in% geneS)] <- ""
########	# same for size. those NOT IN interL set size = 0.5
########	V(connCompRed)$size <- (degGn/10) + 0.5
########	V(connCompRed)$size[!(V(connCompRed)$symbol %in% geneS)] <- 0.5
########	
########	# compute layout
########	set.seed(131)
########	aR <- vcount(connCompRed)^2
########	layO <- layout_with_lgl(connCompRed, maxiter=10000, maxdelta=vcount(connCompRed),  area=aR, coolexp=1.5, repulserad=aR*vcount(connCompRed), cellsize=sqrt(sqrt(aR)), root=NULL)
########	
########	# and print network
########	png(paste0(saveNets, "04_01_elena_selection.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
########		plot(igraph::simplify(connCompRed), vertex.size=V(connCompRed)$size, vertex.label=V(connCompRed)$symbol, vertex.label.cex=0.3, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5], layout=layO)
########	dev.off()
########	
	################################## find communities, from the connected component
	
	# get the biggest connected component from the originalNet
	bG <- which.max(unlist(lapply(decompose(originalNet), vcount))) 
	connComp <- decompose(originalNet)[[bG]]
	
	# compute edge betwenneess, only for community detection
	E(connComp)$eBetw <- edge_betweenness(connComp)

	# compute communities using louvain
	louvain_partition <- igraph::cluster_louvain(connComp, weights=E(connComp)$eBetw)
	connComp$community <- louvain_partition$membership
	names(connComp$community) <- V(connComp)$name
	
	
	# define a function to weight edges in the network, to obtain a nice plot
	# https://stackoverflow.com/questions/16390221/how-to-make-grouped-layout-in-igraph
	
	weightCommunity <- function(row, membership, weigth.within, weight.between){
		# row <- get.edgelist(connComp) ; membership <- connComp$community ; weigth.within <- 1000 ; weight.between <- 1 ;
	
		if(as.numeric(membership[which(names(membership) == row[1])]) == as.numeric(membership[which(names(membership) == row[2])])){
			weight <- weigth.within
		} else {
			weight <- weight.between
		}
		return(weight)
	}
	
	# print network
	set.seed(131)
	# get edge lengths to plot
	E(connComp)$edgesLen <- apply(get.edgelist(connComp), 1, weightCommunity, connComp$community, 5000, 2)
	
	# compute layout considering edges length
	connComp$layout <- layout.fruchterman.reingold(connComp, weights=E(connComp)$edgesLen)
	
	# plot network with communities and edges length
	png(paste0(saveNets, "04_01_communities.png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		plot(louvain_partition, connComp, vertex.label="", vertex.size=3, edge.color="dimgray", edge.width=0.1)
	dev.off()
	
	# search for meaningful enrichments
	database <- org.Mm.eg.db
	#getting ENTREZ IDs to perform the analysis https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
	quiet(listOfGenes <- AnnotationDbi::select(database, keys=V(connComp)$symbol, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL"), all=T)
	#found genes, with ENTREZ IDs
	validID <- listOfGenes[!is.na(listOfGenes$ENTREZID), ]
	
	# enrich the communities using GO bio processes
	geneXcomm <- list()
	goBP <- list()
	commS <- as.character(unique(connComp$community))
	firstOne <- matrix(0, ncol=2, nrow=length(commS))
	for (cmt in seq(1, length(commS), by=1)) {
		geneXcomm[[commS[cmt]]] <- V(connComp)$symbol[which(connComp$community==commS[cmt])]
		gXcEntrez <- validID$ENTREZID[validID$SYMBOL %in% geneXcomm[[commS[cmt]]]]
		# GO bio process
		quiet(goBP[[commS[cmt]]] <- as.data.frame(enrichGO(gXcEntrez, database, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)), all=T)
		#translate entrez in symbol before exporting the table
		symbl <- vector()
		for (gSet in goBP[[commS[cmt]]]$geneID) {
			genes <- unlist(strsplit(gSet, "/"))
			symbl <- append(symbl, paste(validID$SYMBOL[validID$ENTREZ %in% genes], collapse="/"))
		}
		goBP[[commS[cmt]]] <- cbind.data.frame(goBP[[commS[cmt]]], symbol=symbl)
#		write.xlsx(goBP[[commS[cmt]]], paste0(saveNets, "community_", cmt, "_GOBioProc.xlsx"))
		write.table(goBP[[commS[cmt]]], paste0(saveNets, "community_", cmt, "_GOBioProc.csv"), row.names=FALSE, col.names=F, quote=FALSE)
		firstOne[cmt, ] <- c(commS[cmt], goBP[[commS[cmt]]][1, ]$Description)
	}
#}











#########################################	################################## now, find interesting communities and colour accordingly
#########################################	
#########################################	# get original network
#########################################	connComp <- decompose(originalNet)[[1]]
#########################################	
#########################################	connComp$community <- ""
#########################################	connComp$colour <- rep("gray", length(V(connComp)))
#########################################	# find interleukins community
#########################################	for (cnT in seq(1, length(unlist(neighList)), by=1)) {
#########################################		ilP <- which(V(connComp)$name==names(unlist(neighList)[cnT]))
#########################################		connComp$community[ilP] <- "interleukin"
#########################################		connComp$colour[ilP] <- "tomato"
#########################################	}
#########################################	
#########################################	# find ribosomal community
#########################################	ribS <- grep("Rpl|rpl|rps|Rps|Lactb", V(connComp)$symbol)
#########################################	connComp$community[ribS] <- "ribosome" 
#########################################	connComp$colour[ribS] <- "forestgreen"

#########################################	aR <- vcount(interLNet)^2
#########################################	layO <- layout_with_lgl(connComp, maxiter=10000, maxdelta=vcount(connComp),  area=aR, coolexp=1.5, repulserad=aR*vcount(connComp), cellsize=sqrt(sqrt(aR)), root=NULL)
#########################################	
#########################################	plot(igraph::simplify(connComp), vertex.size=2, vertex.label="", vertex.label.cex=0.3, vertex.color=connComp$colour, vertex.frame.color="white", vertex.label.color="black", vertex.color=brewer.pal(12, "Paired")[5], layout=layO)
