#' ---
#' title: "Cleaning data, for AD"
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

# Load libraries
library("bluster")
library("celldex")
library("cluster")
library("compositions")
library("cowplot")
library("DropletUtils")
library("ggplot2")
library("scuttle")
library("gridExtra")
library("limma")
library("org.Mm.eg.db")
library("pheatmap")
library("plyr")
library("RColorBrewer")
library("scater")
library("scran")
library("Seurat")
library("SingleCellExperiment")
library("SingleR")
library("compositions")
library("AnnotationHub")

source("00_functions.R")

# about [setting seeds](https://stackoverflow.com/questions/18696331/how-many-times-to-use-set-seed123)

# IMPORT DATA  

# set paths
threeXDir <- "/home/gabriele/Desktop/nuovo_sequenziamento/run_441/3X_raw/raw_feature_bc_matrix/"
wtDir <- "/home/gabriele/Desktop/nuovo_sequenziamento/run_441/WT_raw/raw_feature_bc_matrix/"

rawPath <- "/home/gabriele/Desktop/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"

saveIn <- paste0(rawPath, "ANALYSIS_FINAL_441/")
saveImg <- paste0(saveIn, "figures/")

#check if saveIn directory exist, otherwise create it
ifelse(dir.exists(saveIn), TRUE, dir.create(saveIn, recursive=T))
ifelse(dir.exists(saveImg), TRUE, dir.create(saveImg, recursive=F))

print(paste0("analysing exp ", expN))

# load data
rawSce <- list()
rawSce[["3xTG"]] <- read10xCounts(threeXDir, col.names = T)
rawSce[["WT"]] <- read10xCounts(wtDir, col.names = T)

# set palette with 12 colours for colouring plots
colorsPal <- brewer.pal(12, "Paired")

#store gene names
featNames <- colnames(rawSce[["3xTG"]])[-1]

#find proper database for mouse
ah <- AnnotationHub()
queryAh <- query(ah, "Mus musculus")
#there are many! choose one, e.g.
#queryAh$title[grep("Ensembl", queryAh$title)]
#queryAh$ah_id[grep("Ensembl", queryAh$title)]

#load database "Ensembl 103 EnsDb for Mus musculus"
musMusculus <- AnnotationHub()[["AH89457"]]

# start labelling the images
lbl <- 1

# remove unwanted genes from rawSce
for (n in names(rawSce)) {

	#' QUALITY CONTROL STARTS BY REMOVING EMPTY DROPLETS

	#' [RETRIEVING EMPTY DROPLETS](http://bioconductor.org/books/release/OSCA/droplet-processing.html)  
	#' cellranger remove empty droplets [by itself](https://support.bioconductor.org/p/123554/#123562)
	#' since we are using filtered data, filtered by using RSEC and DBEC, the data should contain only putative cells. Hence, there is no need to filter out empty droplets  
	#' in any case BD Rhapsody is based on microwell cells capture, which differs from droplet-based cell capture 
	# plot to search for empty droplets: http://bioconductor.org/books/release/OSCA/droplet-processing.html#background-6

	br.out <- barcodeRanks(counts(rawSce[[n]]), lower = 100, fit.bounds = c(1000, 15000))

	png(paste0(saveImg, "01_", lbl, "_barcodeRanksPlot_", n, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(br.out$rank, br.out$total, log="xy", xlab="Rank", ylab="Total")
		o <- order(br.out$rank)
		lines(br.out$rank[o], br.out$fitted[o], col="red")
		abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
		abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
		legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), legend=c("knee", "inflection"))
	dev.off()
	
	lbl <- lbl + 1

	e.out <- emptyDrops(counts(rawSce[[n]]), lower = metadata(br.out)$inflection)

	is.cell <- e.out$FDR <= 0.001
	sum(is.cell, na.rm=TRUE)
	
	write.table(table(Limited=e.out$Limited, Significant=is.cell), "limited_x_significant.csv", quote=F)

	png(paste0(saveImg, "01_", lbl, "_barcodeRanksPlot_", n, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(e.out$Total, -e.out$LogProb, col=ifelse(is.cell, "red", "black"), xlab="Total UMI count", ylab="-Log Probability")
	dev.off()
	lbl <- lbl + 1
	
	# finally remove empty droplets, or keep those cells that are actual cells
	rawSce[[n]] <- rawSce[[n]][, which(e.out$FDR <= 0.001)]
	
	median(apply(counts(rawSce[[n]]), 2, function(x) sum(x>0)))

	# get gene names
	mappedIds <- mapIds(musMusculus, keys=rowData(rawSce[[n]])$ID, keytype="GENEID", column="SYMBOL")
	
	# discard unkwnown genes
	rawSce[[n]] <- rawSce[[n]][which(is.na(mappedIds)==FALSE), ]
	
	rownames(rawSce[[n]]) <- mappedIds[which(is.na(mappedIds)==FALSE)]
	
	# remove RIKEN GENES AS DESCRIBED IN KIPNIS 2020 AND HERE: https://bioinformatics.stackexchange.com/questions/3205/what-are-riken-genes
	
	#' as mentioned [here](https://bioinformatics.stackexchange.com/questions/3197/how-to-filter-ribosomal-rna-from-scrna-seq-data){target="_blank"}
	#' and using [this](http://ribosome.med.miyazaki-u.ac.jp/rpg.cgi?mode=orglist&org=Mus%20musculus&type=m){target="_blank"} as reference  

	# find riken genes
	riks <- grep("Rik", rownames(rawSce[[n]]))
	#which percentage??
	print(paste0("percentage of RIKEN genes in ", n, " ", round((length(riks)/length(rownames(rawSce[[n]])))*100, digits=2), "%"))	

	# remove riken genes, if any
	if (length(riks) > 0) {
		rawSce[[n]] <- rawSce[[n]][-riks, ]
	}
	
	# and ribosomal proteins
	ribs <- grep("Rpl|rpl|rps|Rps|Lactb", rownames(rawSce[[n]]))
	#which percentage??
	print(paste0("percentage of ribosomal genes in ", n, " ", round((length(ribs)/length(rownames(rawSce[[n]])))*100, digits=2), "%"))
	
	# remove ribosomal genes, if any
	if (length(ribs) > 0) {
		rawSce[[n]] <- rawSce[[n]][-ribs, ]
	}
}

#' REMOVE LOW QUALITY CELLS  

#' this is done one experiment at a time as mentioned [here](http://bioconductor.org/books/release/OSCA/quality-control.html#qc-batch){target="_blank"}  

#' More complex studies may involve batches of cells generated with different experimental parameters (e.g., sequencing depth). In such cases, the adaptive strategy should be applied to each batch separately.  

#' computing quality/control metrics using [adaptive thresholds](https://osca.bioconductor.org/quality-control.html#quality-control-outlier){target="_blank"}  

# declare variable to store the purged experiments
cleanSce <- list()
# and dataframes with the info
dfs <- list()

# loop through the clean experiments, to remove below-threshold cells
for (n in names(rawSce)) {

	# find chromosomal genes
	miTc <- grep("mt-", rownames(rawSce[[n]]))
	# set chromosomal location vector
	chrLoc <- rep("", nrow(rawSce[[n]]))
	chrLoc[miTc] <- "MT"
	chrLoc[!miTc] <- "unknown"
	
	# finally add chromosomal locations
	rowData(rawSce[[n]])$chrLoc <- chrLoc

	# temporary variable for sce, so that rawSce do not change after cells removal
	current <- rawSce[[n]]

	# To obtain an adaptive threshold, we assume that most of the dataset consists
	# of high-quality cells. We then identify cells that are outliers for the
	# various QC metrics, based on the median absolute deviation (MAD) from the
	# median value of each metric across all cells.
	
	# set seed
	set.seed(113)
	
	# retrieve mithocondrial genes, that are labelled as "chrM"
	is.mito <- which(rowData(current)$chrLoc=="MT")
	# compute stats and store them as column data
	current <- addPerCellQC(current, subsets=list(Mito=is.mito))
	# find outliers using the adaptive thresholding over the just computed stats
	qc.mito <- isOutlier(current$subsets_Mito_percent, type="higher")
	qc.lib <- isOutlier(current$sum, log=T, type="lower")
	qc.nexprs <- isOutlier(current$detected, log=T, type="lower")
	
	# get info to plot. see below!
	dis.mito <- which(qc.mito==T)
	dis.lib <- which(qc.lib==T)
	dis.nexprs <- which(qc.nexprs==T)
	
	# A cell that is an outlier for any of these metrics is considered to be of low quality and discarded
	# TRUE means discard, FALSE means keep
	discard <- qc.lib | qc.nexprs | qc.mito
	dis.discard <- which(discard==TRUE)
	
	dfs[[n]] <- DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs), MitoProp=sum(qc.mito), Total=sum(discard))
	
	# IMPORTANT NOTE: there is no need to perform the EmptyDroplet retrieval
	# See here: #https://support.bioconductor.org/p/123554/
	
	# IMPORTANT NOTE: WE ARE ABOUT TO GENERATE a subsetted SingleCellExperiment
	# that we would use for downstream analyses
	
	# filter low quality cells and keep those cells we DON'T want to discard: as mentioned above,
	# we must use !discard to keep those that were FALSE
	cleanSce[[n]] <- current[, !discard]
	
	# evaluate discarded cells: http://bioconductor.org/books/release/OSCA/quality-control.html#quality-control-plots
	
	# Another useful diagnostic involves plotting the proportion of mitochondrial counts against some of the other QC metrics.
	# The aim is to confirm that there are no cells with both large total counts and large mitochondrial counts, to ensure that
	# we are not inadvertently removing high-quality cells that happen to be highly metabolically active (e.g., hepatocytes).
	# We demonstrate using data from a larger experiment involving the mouse brain (Zeisel et al. 2015); in this case, we do not
	# observe any points in the top-right corner in Figure 6.3 that might potentially correspond to metabolically active, undamaged cells.
	
	# create plot
	p <- plotColData(current, x="sum", y="subsets_Mito_percent", colour_by=I(discard)) +
		scale_x_log10() +
		theme(plot.title=element_text(hjust=0.5)) +
		ggtitle(paste(n))
	
	png(paste0(saveImg, "01_02_filteredCells.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(p)
	dev.off()
	
	# set seed
	set.seed(113)
	
	# analysing gene expression using the purged single cell experiments
	lost <- calculateAverage(counts(current)[, discard])
	kept <- calculateAverage(counts(current)[, !discard])
	
	logged <- edgeR::cpm(cbind(lost, kept), log=TRUE, prior.count=2)
	logFC <- logged[,1] - logged[,2]
	abundance <- rowMeans(logged)
	
	# If the discarded pool is enriched for a certain cell type, we should observe increased expression of the corresponding marker genes.
	# No systematic upregulation of genes is apparent in the discarded pool in Figure 6.5, suggesting that the QC step did not inadvertently filter out a cell type.
	
	# plotting the discarded cells, coloured by type
	png(paste0(saveImg, "01_03_geneExpression.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		plot(abundance, logFC, xlab="Average count", ylab="Log-FC (lost/kept)", pch=16, main=n, xlim=c(9, 18))
			points(abundance[dis.mito], logFC[dis.mito], col=colorsPal[1], pch=16)
			points(abundance[dis.nexprs], logFC[dis.nexprs], col=colorsPal[4], pch=16)
			points(abundance[dis.lib], logFC[dis.lib], col=colorsPal[6], pch=16)
	dev.off()
	
	# add condition, i.e. WT or 3xTG, to colData
	colData(cleanSce[[n]])$condition <- n

###########################################	# finally, remove doublets
###########################################	madS <- 5
###########################################	cleanSce[[n]] <- remove_high_nGene_libSize(cleanSce[[n]], madS)
}
 
# finally, join the two exps into the same singleCellExp to go on analysing the full dataset
# http://bioconductor.org/books/release/OSCA/integrating-datasets.html#batch-diagnosis
rowData(cleanSce[["3xTG"]]) <- rowData(cleanSce[["WT"]])
# and merge
mergedSce <- cbind(cleanSce[["3xTG"]], cleanSce[["WT"]])

# now, remove genes with zero counts
to.keep <- rowSums(assay(mergedSce, "counts")) > 0
mergedSce <- mergedSce[to.keep, ]

# and rename duplicate gene names and cell names
colnames(mergedSce) <- make.unique(colnames(mergedSce))
rownames(mergedSce) <- make.unique(rownames(mergedSce))

# evaluate combined gene expression
avgCounts <- calculateAverage(mergedSce, size.factors=NULL)
png(paste0(saveImg, "01_04_combined_gene_expression.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
hist(log10(avgCounts), breaks=100, main="combined gene expression", col=colorsPal[5], xlab=expression(Log[10] ~ "average count"))
dev.off()

#' [NORMALISATION](http://bioconductor.org/books/release/OSCA/normalization.html#motivation){target="_blank"}  

#' scaling normalization involves dividing all counts for each cell by a cell-specific scaling factor, often called a “size factor” (Anders and Huber 2010). The assumption here is that any cell-specific bias (e.g., in capture or amplification efficiency) affects all genes equally via scaling of the expected mean count for that cell. The size factor for each cell represents the estimate of the relative bias in that cell, so division of its counts by its size factor should remove that bias. The resulting “normalized expression values” can then be used for downstream analyses such as clustering and dimensionality reduction  

#' 1) [size factor](http://bioconductor.org/books/release/OSCA/normalization.html#library-size-normalization){target="_blank"}  

#' Library size normalization is the simplest strategy for performing scaling normalization. We define the library size as the total sum of counts across all genes for each cell, the expected value of which is assumed to scale with any cell-specific biases. The “library size factor” for each cell is then directly proportional to its library size where the proportionality constant is defined such that the mean size factor across all cells is equal to 1.  

# set seed
set.seed(113)

# TODO:
# should I do this normalisation on mergedSce??
libSF <- librarySizeFactors(mergedSce)

# plot normalised libraries
png(paste0(saveImg, "01_05_norm_libSF_merged.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
hist(log10(libSF), xlab="Log10[Size factor]", col=colorsPal[5], main="Size Factors")
dev.off()

#' 2) [deconvolution](http://bioconductor.org/books/release/OSCA/normalization.html#normalization-by-deconvolution){target="_blank"}  

#' However, single-cell data can be problematic for these bulk normalization methods due to the dominance of low and zero counts. To overcome this, we pool counts from many cells to increase the size of the counts for accurate size factor estimation. Pool-based size factors are then “deconvolved” into cell-based factors for normalization of each cell’s expression profile. This is performed using the calculateSumFactors() function from scran  

# TODO:
# to block, or not to block?
# clust <- quickCluster(mergedSce)

# set seed
set.seed(113)

# Library size normalization with deconvolution and one blocking variable
clust <- quickCluster(mergedSce, block=factor(mergedSce$condition))

#' We use a pre-clustering step with quickCluster() where cells in each cluster are normalized separately and the size factors are rescaled to be comparable across clusters. This avoids the assumption that most genes are non-DE across the entire population - only a non-DE majority is required between pairs of clusters, which is a weaker assumption for highly heterogeneous populations. By default, quickCluster() will use an approximate algorithm for PCA based on methods from the irlba package. The approximation relies on stochastic initialization so we need to set the random seed for reproducibility.  

#' [min.mean](https://rdrr.io/bioc/scran/man/computeSumFactors.html)  

# By default, we set min.mean to 1 for read count data and 0.1 for UMI data. The exact values of these defaults are more-or-less arbitrary and are retained for historical reasons. The lower threshold for UMIs is motivated by (i) their lower count sizes, which would result in the removal of too many genes with a higher threshold; and (ii) the lower variability of UMI counts, which results in a lower frequency of zeroes compared to read count data at the same mean. We use the median library size to detect whether the counts are those of reads (above 100,000) or UMIs (below 50,000) to automatically set min.mean. Mean library sizes in between these two limits will trigger a warning and revert to using min.mean=0.1.

# minmean <- summary(avgCounts)[3]

# set minmean equal to 0.5 since this is the value used in the previous analysis, i.e. the one in the manuscript
minmean <- 0.5

# set seed
set.seed(113)

deconvSF <- calculateSumFactors(mergedSce, cluster=clust, min.mean=minmean)

#' By default, we set min.mean to 1 for read count data and 0.1 for UMI data. The exact values of these defaults are more-or-less arbitrary and are retained for historical reasons. The lower threshold for UMIs is motivated by (i) their lower count sizes, which would result in the removal of too many genes with a higher threshold; and (ii) the lower variability of UMI counts, which results in a lower frequency of zeroes compared to read count data at the same mean. We use the median library size to detect whether the counts are those of reads (above 100,000) or UMIs (below 50,000) to automatically set min.mean. Mean library sizes in between these two limits will trigger a warning and revert to using min.mean=0.1.  

#' If clusters is specified, filtering by min.mean is performed on the per-cluster average during within-cluster normalization, and then on the (library size-adjusted) average of the per-cluster averages during between-cluster normalization.  

# and plot!
png(paste0(saveImg, "01_06_libSF_Deconv.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
plot(libSF, deconvSF, xlab="Library size factor", ylab="Deconvolution size factor", log='xy', pch=16, col=clust)
abline(a=0, b=1, col="red")
dev.off()

#' [apply size factors](http://bioconductor.org/books/release/OSCA/normalization.html#applying-the-size-factors){target="_blank"}  

#' Once we have computed the size factors, we use the logNormCounts() function from scater to compute normalized expression values for each cell. This is done by dividing the count for each gene/spike-in transcript with the appropriate size factor for that cell. The function also log-transforms the normalized values, creating a new assay called "logcounts". These log-values will be the basis of our downstream analyses in the following chapters.  

#' as mentioned [here](http://bioconductor.org/books/release/OSCA/normalization.html#normalization-transformation){target="_blank"} and [here](https://rdrr.io/bioc/scran/man/computeSumFactors.html){target="_blank"}, min.mean is set in order to remove negative-size factors  

#' If too many genes have consistently low counts across all cells, even the pool-based size factors will be zero. This results in zero or negative size factor estimates for many cells. We avoid this by filtering out low-abundance genes using the min.mean argument. This represents a minimum threshold min.mean on the library size-adjusted average counts from calculateAverage.  

# set seed
set.seed(113)

mergedSce <- computeSumFactors(mergedSce, cluster=clust, min.mean=minmean)
mergedSce <- logNormCounts(mergedSce)

# now print sumFactors and library size
png(paste0(saveImg, "01_07_sumFactors_1.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
par(mfrow=c(1, 2))
# sumFactors
plot(mergedSce$total, sizeFactors(mergedSce), log="xy", col=factor(mergedSce$condition), xlab="Library size", ylab="Size factor", main="Colored by sample")
legend("topleft", pch=1, col=1:4, legend=as.character(levels(factor(mergedSce$condition))))
# library size
plot(mergedSce$total, sizeFactors(mergedSce), log="xy", col=rainbow(n=48)[clust], xlab="Library size", ylab="Size factor", main="Colored by cluster")
dev.off()

# and plotting
png(paste0(saveImg, "01_08_librarySize_1.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
hist(log10(as.matrix(logcounts(mergedSce))), xlab="Log10[Size factor]", col=colorsPal[5], main="logNormCounts")
dev.off()

#' [CELL TYPE ANNOTATION](http://bioconductor.org/books/release/OSCA/cell-type-annotation.html#cell-type-annotation){target="_blank"}  

#' Compiled references over the full dataset  

# 2) with [immgen](https://bioconductor.org/packages/release/data/experiment/vignettes/celldex/inst/doc/userguide.html#immunological-genome-project-immgen){target="_blank"} references

immgen <- celldex::ImmGenData()

#' The ImmGen reference consists of microarray profiles of pure mouse immune cells from the project of the same name (Heng et al. 2008). This is currently the most highly resolved immune reference - possibly overwhelmingly so, given the granularity of the fine labels.  

# store the info into mergedSce
mergedSce$immpred_main <- SingleR(test=mergedSce, ref=immgen, labels=immgen$label.main)
mergedSce$immpred_fine <- SingleR(test=mergedSce, ref=immgen, labels=immgen$label.fine)

# and as pruned labels only
mergedSce$pruned_main <- mergedSce$immpred_main$pruned.labels
mergedSce$pruned_fine <- mergedSce$immpred_fine$pruned.labels

# table(immPred$labels)

# finally, store the normalised object with assigned cell type
saveRDS(object=mergedSce, file=paste0(saveIn, "QCed_data_", expN, ".Rds"))
