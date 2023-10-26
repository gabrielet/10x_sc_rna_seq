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

print("")
print(paste0("analysing all the genes for exp ", expN, " clustered with Zinbwave"))
print("")

# load data, all cells
zbWaved <- readRDS(paste0(saveIn, paste0("full_zbwaved_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds")))

# transform SCE back to Seurat to use the clustering algorithms
zbwSeurat <- as.Seurat(zbWaved, data="counts")

# SNN graph
zbwSeurat <- FindNeighbors(object=zbwSeurat, reduction="zinbwave", dims=1:20)

# SLM clustering with different resolution
seurat_clusters_0.1 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.1, random.seed=1)$seurat_clusters
seurat_clusters_0.3 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.3, random.seed=1)$seurat_clusters
seurat_clusters_0.5 <- FindClusters(object=zbwSeurat, algorithm=3, resolution=0.5, random.seed=1)$seurat_clusters

slmClusters <- as.matrix(data.frame(slm_0.1=seurat_clusters_0.1, slm_0.3=seurat_clusters_0.3, slm_0.5=seurat_clusters_0.5))

#store clustering results in zinbSce
zbWaved$slm_0.1 <- as.factor(slmClusters[, 1])
zbWaved$slm_0.3 <- as.factor(slmClusters[, 2])
zbWaved$slm_0.5 <- as.factor(slmClusters[, 3])

# new labels that compress the fine name into something similar to main labels
zbWaved$fine_to_main <- gsub(" ", "", unlist(lapply(strsplit(zbWaved$pruned_fine, "\\("), `[[`, 1)))

# set colours for plotting. remove RED and LIGHT-BLUE from the cols25() palette
palettE <- unname(cols25()[-c(2, 8)][1:(length(unique(zbWaved$pruned_main)))])
# set a nice red for Tgd
names(palettE) <- unique(zbWaved$pruned_main)

# loop through main labels
colR <- vector()

# for each entry in palettE
for (cPos in 1:length(names(palettE))) {
	# get the color for each cell
	colR[grep(names(palettE)[cPos], zbWaved$pruned_main)] <- palettE[cPos]
}

# at this point use the fine labels to find neutrophils and CD8 and colour accordingly
colR[grep("T.8|T.CD8", zbWaved$pruned_fine)] <- "#E41A1C"

# plot 
png(paste0(saveImg, "03_04_zinbwave_tSNE_ALL_CLUSTERS_con_etichette_", expN, "_K_", kZinb, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(zbWaved, "TSNE", colour_by="pruned_main", text_by="slm_0.5", theme_size=3) + scale_color_manual(values=colR) + guides(color=guide_legend("Cell type")) + theme(legend.position="none") + labs(title="Clustered by cell type"))
dev.off()

# plot
png(paste0(saveImg, "03_06_zinbwave_tSNE_ALL_CLUSTERS_", expN, "_K_", kZinb, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotReducedDim(zbWaved, "TSNE", colour_by="pruned_main", text_by="pruned_main") + scale_color_manual(values=colR) + guides(color=guide_legend("Cell type")) + theme(legend.position="none") + labs(title="Clustered by cell type"))
dev.off()

# plot condition only
png(paste0(saveImg, "03_01_zinbwave_tSNE_ALL_CLUSTERS_condition_only_", expN, "_K_", kZinb, ".png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotTSNE(zbWaved, colour_by="condition"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "03_01_slm_0.5_ALL_CLUSTERS_BY_CLUSTER_", expN, "_K_1000.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotTSNE(zbWaved, colour_by="slm_0.5") + theme(legend.position="none") + labs(title="Clusters"))
dev.off()

# plot tSNE by cluster
png(paste0(saveImg, "03_01_slm_0.5_ALL_CLUSTERS_BY_CLUSTER_con_etichette", expN, "_K_1000.png"), type="cairo", units="in", width=8, height=6, pointsize=12, res=300)
	print(scater::plotTSNE(zbWaved, colour_by="slm_0.5", text_by="slm_0.5") + theme(legend.position="none") + labs(title="Clusters"))
dev.off()

##########################################################################################################################################################################################

# finally, export the data
saveRDS(object=zbWaved, file=paste0(saveIn, "slm_clustered_data_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))

# FINALLY, add clusters info to cleanSce
cleanSce <- readRDS(paste0(saveIn, paste0("QCed_data_", expN, ".Rds")))

cleanSce$slm_0.1 <- zbWaved$slm_0.1
cleanSce$slm_0.3 <- zbWaved$slm_0.3
cleanSce$slm_0.5 <- zbWaved$slm_0.5

# then export cleanSce with clustering info
saveRDS(object=cleanSce, file=paste0(saveIn, "QCed_data_with_clusters_", expN, "_K_", kZinb, "_top_", perC, ".Rds"))

