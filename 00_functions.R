#' ---
#' title: "A collection of functions that are used to plot and enrich"
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

# function that remove unwanted cells using reads and genes counts

remove_high_nGene_libSize <- function(sceObj, log = FALSE, type = "higher"){
	libSize <- colSums(counts(sceObj))
	no_keep_libSize <- isOutlier(libSize, nmads=5, log = log, type = type)
	nGenes <- colSums(counts(sceObj) > 0)
	no_keep_nGenes <- isOutlier(nGenes, nmads=5, log = log, type = type)
	return(sceObj[, !no_keep_libSize | !no_keep_nGenes])
}

# function that plots the histograms and heatmap together

plot_membership <- function(df, cond, col) {

#	df <- df; cond <- conD[2]; col <- colorsPal[2];
	
	df <- df[df$Condition == cond, c(1:2, ncol(df))]
	p21 <- ggplot(df,aes(x=Cluster, y=CellType)) +
		geom_tile(width=0.9, height=0.9, aes(alpha=log(Frequency+1)), fill=col) +
		# coord_fixed() +
		theme_minimal() +
		geom_text(aes(label=ifelse(Frequency>0, Frequency, ""), alpha=log(Frequency+1)), color="white") +
		theme(legend.position="none") +
		scale_y_discrete(limits=lims) +
		scale_fill_manual(values=col)
	df_tot_CellType <- ddply(df, .variables=~ CellType, .fun=function(x) sum(x$Frequency))
	myT <- max(df_tot_CellType$V1)
	colnames(df_tot_CellType) <- c("CellType", "Frequency")
	# computing values as percentages, to keep the histograms away from too big values (that may happen)
	df_tot_CellType$Frequency <- round(df_tot_CellType$Frequency/sum(df_tot_CellType$Frequency) * 100, digits = 1)
	p22 <- ggplot(df_tot_CellType,aes(x=CellType, y=Frequency)) +
		geom_col(aes(alpha=log(Frequency+1)), fill=col) +
		theme_minimal() +
		coord_flip() +
		geom_text(aes(y=Frequency + 30, label=ifelse(Frequency>0, Frequency, ""), alpha=log(Frequency+1)), hjust="left") +
		theme(legend.position="none",
			axis.title=element_blank(),
			axis.text=element_blank(), panel.grid=element_blank()) +
		scale_x_discrete(limits=lims) +
		scale_y_continuous(limits=c(0, myT)) +
		scale_fill_manual(values=col)
	df_tot_Cluster <- ddply(df, .variables=~ Cluster, .fun=function(x) sum(x$Frequency))
	myC <- max(df_tot_Cluster$V1)
	colnames(df_tot_Cluster) <- c("Cluster","Frequency")
	# computing values as percentages, to keep the histograms away from too big values (that may happen)
	df_tot_Cluster$Frequency <- round(df_tot_Cluster$Frequency/sum(df_tot_Cluster$Frequency) * 100, digits = 1)
	p11 <- ggplot(df_tot_Cluster,aes(x=Cluster, y=Frequency)) +
		geom_col(aes(alpha=log(Frequency+1)), fill=col) +
		theme_minimal() +
		geom_text(aes(y=Frequency + 150, label=ifelse(Frequency>0, Frequency, ""), alpha=log(Frequency+1))) +
		theme(legend.position="none",
			axis.title=element_blank(),
			axis.text=element_blank(), panel.grid=element_blank()) +
		scale_y_continuous(limits=c(0, myC)) +
		scale_fill_manual(values=col)
	empty <- ggplot() +
		geom_point(aes(1,1), colour="white") +
		theme(axis.ticks=element_blank(),
			panel.background=element_blank(),
			axis.text.x=element_blank(),
			axis.text.y=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank())
	plot_grid(p11, empty, p21, p22, rel_heights=c(0.2,1), rel_widths=c(1,0.4), ncol=2,nrow=2, align="hv",axis="lrtb")
}

# plotting upset plots

# using this version: https://github.com/hms-dbmi/UpSetR/issues/85 of fromList() function

fromList <- function (input) {
	# Same as original fromList()...
	elements <- unique(unlist(input))
	data <- unlist(lapply(input, function(x) {x <- as.vector(match(elements, x))}))
	data[is.na(data)] <- as.integer(0)
	data[data != 0] <- as.integer(1)
	data <- data.frame(matrix(data, ncol = length(input), byrow = F))
	data <- data[which(rowSums(data) != 0), ]
	names(data) <- names(input)
	# ...Except now it conserves its original value names!
	row.names(data) <- elements
	# finally, return results
	return(data)
}

# analysing clusters

################## OPTION ONE FOR NOT-BATCHED AD DATA

optionOne <- function(subSet, exP, saveGReg, myT, cnD, cTy, thresHD) {

#	subSet <- cTSubSet; exP <- expN; saveGReg <- saveUpDown; myT <- personalT; cnD <- conD; cTy <- cellLab[cT]; thresHD <- thresholD

	cTSce <- list()

	# then get WT cellType
	cTSce[[cnD[1]]] <- subSet[, which(subSet$condition==cnD[1])]
	# and 3xTG
	cTSce[[cnD[2]]] <- subSet[, which(subSet$condition==cnD[2])]

	#print info about number of cells for each condition
	print(paste0("in cluster ", cS, " there are ", ncol(cTSce[[cnD[1]]]), " ", cTy, " for ", cnD[1]))
	print(paste0("and ", ncol(cTSce[[cnD[2]]]), " ", cTy, " for ", cnD[2]))

	# this if else check is useful here and not for the two other options since here all the cellType cells are considered
	if (ncol(cTSce[[cnD[1]]]) > myT && ncol(cTSce[[cnD[2]]]) > myT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")

		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cnD[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering cnD[1] as the chosen condition
		cC <- cnD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]		

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells in cluster ", cS))
		print(paste0("for the comparison of ", cC, " over ", cnD[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]

		write.table(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", exP, "_UP_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", exP, "_UP_REG.xlsx"))

		write.table(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", exP, "_DOWN_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", exP, "_DOWN_REG.xlsx"))

		write.table(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", exP, "_up_and_down.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_cluster_", cS, "_exp_", exP, "_up_and_down.xlsx"))

		# store the genes for this cluster in the list
		optOneGenes <- finalGenes$gName

	} else {
		# set return to NULL
		optOneGenes <- NULL
		print(paste0("not enough cells to compare the two conditions for ", cTy, " in cluster ", cS))
	}
	# finally, return the results
	return(optOneGenes)
}

################## OPTION TWO FOR NOT-BATCHED AD DATA

optionTwo <- function(subSet, myT, cnD, cTy, thresHD) {

	cTSce <- list()

	# then get WT cellType
	cTSce[[cnD[1]]] <- subSet[, which(subSet$condition==cnD[1])]
	# and 3xTG
	cTSce[[cnD[2]]] <- subSet[, which(subSet$condition==cnD[2])]

	#print info about number of cells for each condition
	print(paste0("in cluster ", cS, " there are ", ncol(cTSce[[cnD[1]]]), " cells for ", cnD[1]))
	print(paste0("and ", ncol(cTSce[[cnD[2]]]), " cells for ", cnD[2]))

	if (ncol(cTSce[[cnD[1]]]) > myT && ncol(cTSce[[cnD[2]]]) > myT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")

		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cnD[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering cnD[1] as the chosen condition
		cC <- cnD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for cluster ", cS))
		print(paste0("for the comparison of ", cC, " over ", cnD[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]
		
		# store the genes for this cluster in the list
		optTwoGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optTwoGenes <- NULL
		print(paste0("not enough cells to compare the two conditions for cluster ", cS))
	}
	# finally, return the results
	return(optTwoGenes)
}

################## OPTION THREE FOR NOT-BATCHED AD DATA

optionThree <- function(subSet, exP, saveGReg, myT, cnD, cTy, thresHD) {

	cTSce <- list()

	# then get WT cellType
	cTSce[[cnD[1]]] <- subSet[, which(subSet$condition==cnD[1])]
	# and 3xTG
	cTSce[[cnD[2]]] <- subSet[, which(subSet$condition==cnD[2])]

	#print info about number of cells for each condition
	print(paste0("in the full dataset there are ", ncol(cTSce[[cnD[1]]]), " ", cTy, " for ", cnD[1]))
	print(paste0("and ", ncol(cTSce[[cnD[2]]]), " ", cTy, " for ", cnD[2]))

	if (ncol(cTSce[[cnD[1]]]) > myT && ncol(cTSce[[cnD[2]]]) > myT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")
		
		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering cnD[1] as the chosen condition
		cC <- cnD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells"))
		print(paste0("for the comparison of ", cC, " over ", cnD[2]))

		# and export data
		finalGenes <- cbind.data.frame(gName=rownames(selfMade), selfMade)

		# order by fold change, which is the 6th column
		finalOrdered <- list()
		finalOrdered[["UpAndDown"]] <-  finalGenes[order(finalGenes[, 6], decreasing=T), ]

		# get up and down, separated
		finalOrdered[["Up"]] <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]>=0), ]
		dTmp <- finalOrdered[["UpAndDown"]][which(finalOrdered[["UpAndDown"]][, 6]<0), ]
		finalOrdered[["Down"]] <- dTmp[order(dTmp[, 6], decreasing=F), ]

		write.table(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_UP_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Up"]], paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_UP_REG.xlsx"))

		write.table(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_DOWN_REG.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["Down"]], paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_DOWN_REG.xlsx"))

		write.table(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_up_and_down.csv"), quote=F, row.names=F, col.names=T, sep="\t")
		write.xlsx(finalOrdered[["UpAndDown"]], paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_up_and_down.xlsx"))
		
		# store the genes for this cluster in the list
		optThreeGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optThreeGenes <- NULL
		print(paste0("not enough ", cTy, " cells to compare the two conditions"))
	}
	return(optThreeGenes)
}

################## OPTION THREE FOR random data

optionThreeRandomSampling <- function(subSet, exP, myT, cnD, cTy, thresHD) {

	cTSce <- list()

	# then get WT cellType
	cTSce[[cnD[1]]] <- subSet[, which(subSet$condition==cnD[1])]
	# and 3xTG
	cTSce[[cnD[2]]] <- subSet[, which(subSet$condition==cnD[2])]

	#print info about number of cells for each condition
	print(paste0("in the full dataset there are ", ncol(cTSce[[cnD[1]]]), " ", cTy, " for ", cnD[1]))
	print(paste0("and ", ncol(cTSce[[cnD[2]]]), " ", cTy, " for ", cnD[2]))

	if (ncol(cTSce[[cnD[1]]]) > myT && ncol(cTSce[[cnD[2]]]) > myT) {

		# then search for marker genes using pairwise by condition so that WT will be compared to 3x (and vice versa)
		# https://rdrr.io/bioc/scran/man/findMarkers.html
		# it is important to note that findMarkers uses logcounts as the default assay, to compute the differentially expressed genes.
		# however when we apply zinbwave, its normalised values are put into the logcounts assay
		# also, only 1000 highly variable genes are used
		
		markersG <- findMarkers(subSet, groups=subSet$condition, pval.type="any", direction="any")
		
		# http://bioconductor.org/books/release/OSCA/marker-detection.html#looking-for-any-differences
		# The default philosophy of findMarkers() is to identify a combination of marker genes that - together - uniquely define one cluster against the rest. 
		# To this end, we collect the top DE genes from each pairwise comparison involving a particular cluster to assemble a set of candidate markers for that cluster.
		# We will demonstrate on cluster 7; the relevant DataFrame contains log2-fold changes of expression in cluster 7 over each other cluster, along with several 
		# statistics obtained by combining p-values across the pairwise comparisons involving 7.
		# chosen <- "7"
		# interesting <- markers.pbmc[[chosen]]

		# hence, if cond[1] stands for "3xTG", we are obtaining the fold change for the comparison of 3xTG over WT

		# filter on false discovery rate considering cnD[1] as the chosen condition
		cC <- cnD[1]
		chosen <- markersG[[cC]]
		selfMade <- chosen[chosen$FDR < thresHD, ]

		print(paste0(nrow(selfMade), " statistical significant genes were found for all the ", cTy, " cells"))
		print(paste0("for the comparison of ", cC, " over ", cnD[2]))

		# store the genes for this cluster in the list
		optThreeGenes <- finalGenes$gName
	} else {
		# set return to NULL
		optThreeGenes <- NULL
		print(paste0("not enough ", cTy, " cells to compare the two conditions"))
	}
	return(optThreeGenes)
}

################## findMarkers with Seurat

fMSeurat <- function(srtSubSet, resList, thD, saveGReg, exP, cTy) {

#	srtSubSet <- seuratSubSet; resList <- seuratMarkers; thD <- thresholD; saveGReg <- saveUpDown; exP <- expN; cTy <- cellLab[cT]

	methods <- c("wilcox", "bimod", "roc", "t", "poisson", "negbinom", "LR", "MAST", "DESeq2")
	
	for (mth in methods) {
		reS <- FindMarkers(srtSubSet, ident.1=conD[[1]], ident.2=conD[[2]], test.use=mth)
		resList[[mth]] <- reS[which(reS$p_val_adj <  thD), ]
		# to export
		toExp <- data.frame(gName=rownames(resList[[mth]]), resList[[mth]])
		write.table(toExp, paste0(saveGReg, cTy, "_ONLY_exp_", exP, "_up_and_down_with_", mth, ".csv"), quote=F, row.names=T, col.names=T, sep="\t")
	}
	return(resList)
}

# perform enrichment using DESeq2
rankMyData <- function(subSet) {

	#subSet <- cTSubSet;
	
	# obtain a list of DEGs
	subSet$condition <- as.factor(subSet$condition)
	quiet(dds <- DESeqDataSetFromMatrix(as.matrix(assays(subSet)$counts), colData=colData(subSet), design=~condition), all=T)
	## Pre filter

	# We perform a minimal pre-filtering to keep only rows that have at least 10 reads total. Note that more strict filtering to increase power is automatically applied via independent filtering on the mean of normalized counts within the results function

	keep <- rowSums(counts(dds)) >= 10
	dds <- dds[keep, ]

	# After this first filter, we remove genes without an EntrezID
	rowData(dds)$Symbol <- rownames(dds)
	quiet(rowData(dds)$Entrezid <- mapIds(x=org.Mm.eg.db, keys=rowData(dds)$Symbol, column="ENTREZID", keytype="SYMBOL"), all=T)

	keep <- !is.na(rowData(dds)$Entrezid)
	dds <- dds[keep, ]

	# set "WT" or "CTRL" as the reference level
	if("WT" %in% dds$condition){
		dds$condition <- factor(dds$condition, levels=c("WT", "3xTG"))
	# else EAE
	} else if ("EAE" %in% dds$condition) {
		dds$condition <- factor(dds$condition, levels=c("CTRL", "EAE"))
	}
	
	## Differential Expression

	# We launch DESeq2
	quiet(ddsDesequed <- DESeq(dds), all=T)
	res <- results(ddsDesequed)

	# A positive value of Log-Fold Change indicates higher expression values in 3X_Alz compared to WT

	# The analysis is performed by:
	# - ranking all genes in the data set;
	# - identifying the rank positions of all members of the gene set in the ranked data set;
	# - calculating an enrichment score (ES) that represents the difference between the observed rankings and that which would be expected assuming a random rank distribution.

	## Create ranks
	
	if (all(rownames(res) == rowData(ddsDesequed)$Symbol)) {
		# Rank all genes based on their statistics
		rkVals <- res$stat
		# and name each gene
		names(rkVals) <- rowData(ddsDesequed)$Entrezid
	} else {
		print("problem with names")
	}
	
	# finally, return results
	return(rkVals)
}

enrichMyData <- function(oData, saveCProf, id, thrS, exP) {

#	oData <- originalData; saveCProf <- cProfPathsList; id <- ID; thrS <- threshS;

	#call the database from which we will get the info
	database <- org.Mm.eg.db
	
	#create a string array with the names of the genes
	geneNames <- as.character(oData$gName)

	#getting ENTREZ IDs to perform the analysis https://davetang.org/muse/2013/12/16/bioconductor-annotation-packages/
	quiet(listOfGenes <- AnnotationDbi::select(database, keys=geneNames, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL"), all=T)

	#found genes, with ENTREZ IDs
	validID <- listOfGenes[!is.na(listOfGenes$ENTREZID), ]

	# GO bio process
	quiet(goBP <- as.data.frame(enrichGO(validID$ENTREZID, database, ont="BP", pAdjustMethod="BH", pvalueCutoff=thrS[1], qvalueCutoff=thrS[2])), all=T)
	#translate entrez in symbol before exporting the table
	symbl <- vector()
	for (gSet in goBP$geneID) {
		genes <- unlist(strsplit(gSet, "/"))
		symbl <- append(symbl, paste(validID$SYMBOL[validID$ENTREZ %in% genes], collapse="/"))
	}
	goBP <- cbind.data.frame(goBP, symbol=symbl)
	# then export the table using the first slot of saveCProf paths list
#	write.table(goBP, paste0(saveCProf[1], "cProf_", id, "_GO_bio_process_", exP, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
	write.xlsx(goBP, paste0(saveCProf[1], "cProf_", id, "_GO_bio_process_", exP, ".xlsx"))
	
	# GO mol function
	quiet(goMF <- as.data.frame(enrichGO(validID$ENTREZID, database, ont="MF", pAdjustMethod="BH", pvalueCutoff=thrS[1], qvalueCutoff=thrS[2])), all=T)		
	
	#translate entrez in symbol before exporting the table
	symbl <- vector()
	for (gSet in goMF$geneID) {
		genes <- unlist(strsplit(gSet, "/"))
		symbl <- append(symbl, paste(validID$SYMBOL[validID$ENTREZ %in% genes], collapse="/"))
	}
	goMF <- cbind.data.frame(goMF, symbol=symbl)
	# then export the table using the second slot of saveCProf paths list
#	write.table(goMF, paste0(saveCProf[2], "cProf_", id, "_GO_mol_function_", exP, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
	write.xlsx(goMF, paste0(saveCProf[2], "cProf_", id, "_GO_mol_function_", exP, ".xlsx"))
	
	# REACTOME pathways
	quiet(reactome <- as.data.frame(enrichPathway(validID$ENTREZID, organism="mouse", pAdjustMethod="BH", pvalueCutoff=thrS[1], qvalueCutoff=thrS[2])), all=T)
	#translate entrez into symbol befor exporting the table
	symbl <- vector()
	for (gSet in reactome$geneID) {
		genes <- unlist(strsplit(gSet, "/"))
		symbl <- append(symbl, paste(validID$SYMBOL[validID$ENTREZ %in% genes], collapse="/"))
	}
	reactome <- cbind.data.frame(reactome, symbol=symbl)
	# then export the table using the third slot of saveCProf paths list
#	write.table(reactome, paste0(saveCProf[3], "cProf_", id, "_REACTOME_", exP, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
	write.xlsx(reactome, paste0(saveCProf[3], "cProf_", id, "_REACTOME_", exP, ".xlsx"))
	
	# KEGG pathways
	quiet(kegg <- as.data.frame(enrichKEGG(validID$ENTREZID, organism="mouse", pAdjustMethod="BH", pvalueCutoff=thrS[1], qvalueCutoff=thrS[2])), all=T)
	#translate entrez into symbol befor exporting the table
	symbl <- vector()
	for (gSet in kegg$geneID) {
		genes <- unlist(strsplit(gSet, "/"))
		symbl <- append(symbl, paste(validID$SYMBOL[validID$ENTREZ %in% genes], collapse="/"))
	}
	kegg <- cbind.data.frame(kegg, symbol=symbl)
	# then export the table using the third slot of saveCProf paths list
#	write.table(kegg, paste0(saveCProf[4], "cProf_", id, "_KEGG_", exP, ".csv"), row.names=F, col.names=T, quote=F, sep="\t")
	write.xlsx(kegg, paste0(saveCProf[4], "cProf_", id, "_KEGG_", exP, ".xlsx"))

	# finally, return results
	return(TRUE)
}

# interestingly, it is always possibe to subset any EnrichmentObject using the information available here:

# https://github.com/YuLab-SMU/enrichplot/issues/17
# http://yulab-smu.top/clusterProfiler-book/chapter13.html

# subsetting is surely useful to show only those terms that are actually interesting to the process under investigation

# some GO
goPlotting <- function(rkDataF, saveGOImg, saveCProf, lbl, fN, cS, exP) {

	# rkDataF <- originalData; saveGOImg <- saveGO; lbl <- lbl;

	gene <- rkDataF$gene
	
	# generate Rankings data structure
	ranks <- rkDataF$deseqStat
	names(ranks) <- rkDataF$gene

	#getting SYMBOLS IDs from ENTREZ
	quiet(gSymb <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(gene), columns=c("SYMBOL", "ENTREZID"), keytype="ENTREZID"), all=T)

	eGO <- enrichGO(gene=gene, OrgDb=org.Mm.eg.db, ont="BP", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05)
	
	# entrez to symbols
	symbs <- vector()
	expGO <- as.data.frame(eGO)
	
	if (nrow(expGO) > 0) {
		for (nr in seq(1, nrow(expGO), by=1)) {
			# get symbols
			symbs <- append(symbs, paste(gSymb$SYMBOL[gSymb$ENTREZID %in% strsplit(expGO[nr, 8], "/")[[1]]], collapse="/"))
		}
		# get final table ready
		expGO <- cbind.data.frame(as.data.frame(eGO), symbol=symbs)
		write.table(expGO, paste0(saveCProf, "cProf_DESEQ_on_rankings_GO_Bio_Process_", exP, "_cluster_", cS, ".csv", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(expGO, paste0(saveCProf, "cProf_DESEQ_on_rankings_GO_Bio_Process_", exP, "_cluster_", cS, ".xlsx"))
	} else {
		print("no enrichment")
	}

	if (nrow(eGO) > 0) {
		lbl <- lbl + 1
		p <- goplot(eGO)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "goplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=15, pointsize=12, res=300)
		plot(p)
		dev.off()

		lbl <- lbl + 1
		p <- dotplot(eGO)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "dotplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		plot(p)
		dev.off()

		# setReadable translate entrezID into symbols
		eGOx <- setReadable(eGO, 'org.Mm.eg.db', 'ENTREZID')
		lbl <- lbl + 1
		p <- cnetplot(eGOx, foldChange=ranks, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "cnetplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=32, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		lbl <- lbl + 1
		p <- heatplot(eGOx, foldChange=ranks)
		png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "heatplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in",width=40, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		# to plot a pairwise plot, there should be at least a couple terms
		if (nrow(eGOx) >= 2) {
			eGO <- pairwise_termsim(eGOx)
			lbl <- lbl + 1
			p <- emapplot(eGO, cex_category=1.5,layout="kk")
			png(paste0(saveGOImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "pairwiseplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=18, pointsize=12, res=300)
			plot(p)
			dev.off()
		}
	} else {
		print("no enrichment")
	}
	# return label
	return(lbl)
}

# some REACTOME
reactomePlotting <- function(rkDataF, saveRctImg, saveCProf, lbl, fN, cS, exP) {
	
	# rkDataF <- originalData; saveRctImg <- saveReactome; lbl <- lbl;
	
	# generate Rankings data structure
	ranks <- rkDataF$deseqStat
	names(ranks) <- rkDataF$gene

	#getting SYMBOLS IDs from ENTREZ
	quiet(gSymb <- AnnotationDbi::select(org.Mm.eg.db, keys=as.character(rkDataF$gene), columns=c("SYMBOL", "ENTREZID"), keytype="ENTREZID"), all=T)

	eRCTM <- enrichPathway(gene=rkDataF$gene, organism="mouse", pAdjustMethod="BH", pvalueCutoff=0.01, qvalueCutoff=0.05, readable=TRUE)
	
	# entrez to symbols
	symbs <- vector()
	expRCT <- as.data.frame(eRCTM)
	
	if (nrow(expRCT) > 0) {
		for (nr in seq(1, nrow(expRCT), by=1)) {
			# get symbols
			symbs <- append(symbs, paste(gSymb$SYMBOL[gSymb$ENTREZID %in% strsplit(expRCT[nr, 8], "/")[[1]]], collapse="/"))
		}
		# get final table ready
		expRCT <- cbind.data.frame(as.data.frame(eRCTM), symbol=symbs)
		write.table(expRCT, paste0(saveCProf, "cProf_DESEQ_on_rankings_Reactome_", exP, "_cluster_", cS, ".csv", sep=""), row.names=F, col.names=T, quote=F, sep="\t")
		write.xlsx(expRCT, paste0(saveCProf, "cProf_DESEQ_on_rankings_Reactome_", exP, "_cluster_", cS, ".xlsx"))
	} else {
		print("no enrichment")
	}

	if (nrow(eRCTM) > 0) {
		lbl <- lbl + 1
		p <- dotplot(eRCTM)
		png(paste0(saveRctImg, "04_0", strsplit(gsub("rankings_genes", "dotplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		plot(p)
		dev.off()
		
		# setReadable translate entrezID into symbols
		eRCTx <- setReadable(eRCTM, 'org.Mm.eg.db', 'ENTREZID')
		lbl <- lbl + 1
		p <- cnetplot(eRCTx, foldChange=ranks, categorySize="pvalue", circular = TRUE, colorEdge = TRUE)
		png(paste0(saveRctImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "cnetplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=32, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()
		
		lbl <- lbl + 1
		p <- heatplot(eRCTx, foldChange=ranks)
		png(paste0(saveRctImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "heatplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in",width=40, height=32, pointsize=12, res=300)
		plot(p)
		dev.off()

		# to plot a pairwise plot, there should be at least a couple pathways
		if (nrow(eRCTx) >= 2) {
			eRCT <- pairwise_termsim(eRCTx)
			lbl <- lbl + 1
			p <- emapplot(eRCT, cex_category=1.5,layout="kk")
			png(paste0(saveRctImg, "04_0", lbl, "_", strsplit(gsub("rankings_genes", "pairwiseplot", fN), "\\.")[[1]][1], ".png"), type="cairo", units="in", width=20, height=18, pointsize=12, res=300)
			plot(p)
			dev.off()
		}
	} else {
		print("no enrichment")
	}
	# return label
	return(lbl)
}
