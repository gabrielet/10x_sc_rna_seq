#' ---
#' title: "SLM clustering for zinbwaved data"
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
library("RColorBrewer")
library("scater")
library("scran")
library("SingleCellExperiment")
library("Seurat")
library("pals")

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

# number of highly variable genes used for zinbwaving the data
perC <- 1000
kZinb <- 20

# select cell types
cellType <- c("T.8|T.CD8", "Neutrophils", "Tgd")#, "Monocytes", "Macrophages", "Microglia", "B cells", "ILC", "NK cells", "NKT", "Mast")
cellLab <-  c("CD8", "Neutrophils", "Tgd")#, "Monocytes", "Macrophages", "Microglia", "B_cells", "ILC", "NK_cells", "NKT", "Mast")

# loop through cell type
for (cT in seq(1, length(cellType), by=1)){

	# load data, all cells
	zbWaved <- readRDS(paste0(saveIn, "full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))
	
	print("")
	print(paste0("analysing all the genes ", cellLab[cT], " for exp ", expN, " clustered with Zinbwave"))
	print("")

	###################################################################### get a specific cell type to perform subclustering ######################################################################
	zbWaved <- zbWaved[, grep(cellType[cT], zbWaved$pruned_fine)]
	###############################################################################################################################################################################################

	# transform SCE back to Seurat to use the clustering algorithms
	zbwSeurat <- as.Seurat(zbWaved, data="counts")

	# SNN graph
	zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="zinbwave", dims=1:20)

	# SLM clustering with different resolution	
	seurat_clusters_0.4 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.4, random.seed=1)$seurat_clusters
	seurat_clusters_0.3 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.3, random.seed=1)$seurat_clusters

	slmClusters <- as.matrix(data.frame(slm_0.4=seurat_clusters_0.4, slm_0.3=seurat_clusters_0.3))

	#store clustering results in zinbSce
	zbWaved$slm_0.4<- as.factor(slmClusters[, 1])
	zbWaved$slm_0.3 <- as.factor(slmClusters[, 2])

	# new labels that compress the fine name into something similar to main labels
	zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

	# plot condition only
	png(paste0(saveImg, "03_subcl_02_zinbwave_tSNE_SUBCLUSTERED_condition_only_", expN, "_K_", kZinb, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(zbWaved, colour_by="condition"))
	dev.off()
	
	# plot clusters only
	png(paste0(saveImg, "03_subcl_02_zinbwave_tSNE_SUBCLUSTERED_by_clusters_", expN, "_K_", kZinb, "_", cellLab[cT], ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
		print(scater::plotTSNE(zbWaved, colour_by="slm_0.3"))
	dev.off()

	# plot
	png(paste0(saveImg, "03_subcl_05_zinbwave_tSNE_SUBCLUSTERED_", expN, "_K_", kZinb, "_", cellLab[cT], ".png"), type="cairo", units="in", width=12, height=9, pointsize=12, res=300)
		print(scater::plotReducedDim(zbWaved, dimred="TSNE", colour_by="pruned_fine") + theme(legend.position="bottom") + labs(title="By cell type") + guides(color=guide_legend("Cell type")))
	dev.off()

	# finally, export the data
	saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, "_K_", kZinb, "_top_", perC, "_", cellLab[cT], "_ONLY.Rds"))

	# FINALLY, add clusters info to cleanSce
	cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

	################################################################################## get celltype only to perform subclustering #####################################################################
	cleanSce <- cleanSce[, grep(cellType[cT], cleanSce$pruned_fine)]
	###################################################################################################################################################################################################

	cleanSce$slm_0.4 <- zbWaved$slm_0.4
	cleanSce$slm_0.3 <- zbWaved$slm_0.3

	# then export cleanSce with clustering info
	saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, "_", cellLab[cT], "_ONLY.Rds"))
}
