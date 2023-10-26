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
rawPath <- "/home/gabriele/Desktop/nuovo_sequenziamento/new_analysis_with_RIBOSOMAL/441/" ; expN <- "441"
 
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
zbWaved <- readRDS(paste0(saveIn, paste0("full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds")))

####################################################################################### get CD8 to perform subclustering #######################################################################################
zbWaved <- zbWaved[, grep("T.8|T.CD8", zbWaved$pruned_fine)]
################################################################################################################################################################################################################

# transform SCE back to Seurat to use the clustering algorithms
zbwSeurat <- as.Seurat(zbWaved, data="counts")

# SNN graph
zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="zinbwave", dims=1:20)

# SLM clustering with different resolution
seurat_clusters_0.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.5, random.seed=1)$seurat_clusters
seurat_clusters_1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=1, random.seed=1)$seurat_clusters
seurat_clusters_1.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=1.5, random.seed=1)$seurat_clusters

slmClusters <- as.matrix(data.frame(slm_0.5=seurat_clusters_0.5, slm_1=seurat_clusters_1, slm_1.5=seurat_clusters_1.5))

#store clustering results in zinbSce
zbWaved$slm_0.5 <- as.factor(slmClusters[, 1])
zbWaved$slm_1 <- as.factor(slmClusters[, 2])
zbWaved$slm_1.5 <- as.factor(slmClusters[, 3])

# new labels that compress the fine name into something similar to main labels
zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

# a summary of clusters content
tbl <- table(zbWaved$slm_1, zbWaved$fine_to_main)

# plot tSNE by condition
png(paste0(saveImg, "03_01_slm_0.5_BY_CONDITION_", expN, "_K_1000_CD8_ONLY.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	scater::plotTSNE(zbWaved, colour_by="condition") + theme(legend.position="right") + labs(title="Conditions")
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "03_01_slm_0.5_BY_CLUSTER_", expN, "_K_1000_CD8_ONLY.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	scater::plotTSNE(zbWaved, colour_by="slm_0.5") + theme(legend.position="none") + labs(title="Clusters")
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "03_01_slm_0.5_BY_CLUSTER_con_etichette", expN, "_K_1000_CD8_ONLY.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	scater::plotTSNE(zbWaved, colour_by="slm_0.5", text_by="slm_0.5") + theme(legend.position="none") + labs(title="Clusters")
dev.off()

# and by cells
png(paste0(saveImg, "03_01_slm_0.5_BY_CELLS_", expN, "_K_1000_CD8_ONLY.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	scater::plotReducedDim(zbWaved, "TSNE", colour_by="pruned_fine") + theme(legend.position="bottom") + labs(title="Cell types") + guides(color=guide_legend("Cell type"))
dev.off()

####################################################################################### commented for subclustering #######################################################################################

###################################### set colours for plotting. remove RED and LIGHT-BLUE from the cols25() palette
#####################################palettE <- unname(cols25()[-c(2, 8)][1:(length(unique(zbWaved$fine_to_main))-1)])
###################################### set a nice red for Tgd
#####################################names(palettE) <- unique(zbWaved$fine_to_main)[!is.na(unique(zbWaved$fine_to_main))]
#####################################tgd <- grep("Tgd", names(palettE))
#####################################palettE[tgd] <- "#E41A1C"
###################################### and a nice light-blue for neutrophils
#####################################pmn <- grep("Neutrophils", names(palettE))
#####################################palettE[pmn] <- "#a6cee3"

###########################################################################################################################################################################################################

#png(paste0(saveImg, "03_01_slm_0.5_clustering_only_main_", expN, "_K_1000.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
#	scater::plotReducedDim(zbWaved, "TSNE", colour_by="fine_to_main") + scale_color_manual(values=palettE) + guides(color=guide_legend("Cell type")) + theme(legend.position="right") + labs(title="Clustered by cell type")
#dev.off()

# finally, export the data
saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, "_K_", kZinb, "_top_", perC, "_CD8_ONLY.Rds"))

# FINALLY, add clusters info to cleanSce
cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

####################################################################################### get CD8 to perform subclustering #######################################################################################
cleanSce <- cleanSce[, grep("T.8|T.CD8", cleanSce$pruned_fine)]
################################################################################################################################################################################################################

cleanSce$slm_0.5 <- zbWaved$slm_0.5
cleanSce$slm_1 <- zbWaved$slm_1
cleanSce$slm_1.5 <- zbWaved$slm_1.5

# then export cleanSce with clustering info
saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, "_CD8_ONLY.Rds"))
