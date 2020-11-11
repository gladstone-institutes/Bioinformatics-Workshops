##remove all data: start from scratch
rm(list = ls())
#Load the libraries.
library(Seurat)
library(muscat)
library(SummarizedExperiment)
library(dplyr)
library(magrittr)
library(purrr)
library(tibble)
library(gdata)
library(tidyr)
source("pbDS_update.R")

raw_data <- read.csv("rawCounts.csv", header = T)
pheno_data <- read.csv("sub_pheno_data.csv", header = T)
print(dim(raw_data))
print(dim(pheno_data))
head(pheno_data)

##randomly introduce batch information for illustrating batch correction procedure
##Individuals 1197 and 1235 are assigned to batch 1
##Individuals 1200 and 1242 are assigned to batch 2
batch <- rep("batch1", nrow(pheno_data))
batch[pheno_data$Individual==1200 | pheno_data$Individual==1242] <- "batch2"
pheno_data <- data.frame(pheno_data, batch)

mm10_genes <- read.csv("mm10_genes.tsv", header=FALSE, sep='\t', stringsAsFactors=FALSE,
                       col.names=c("ensembl_id", "gene_symbol"))

gene_ids <- as.character(raw_data$Geneid)

raw_data <- raw_data[,-1]
row.names(raw_data) <- gene_ids

# Map ENSEMBL Ids to their gene symbols
TempIndices <- match(gene_ids, mm10_genes$ensembl_id)
raw_data <- raw_data[!is.na(TempIndices), ]
CheckIds <- row.names(raw_data)[1:5]
NonUniqueGeneSymbols <-  mm10_genes$gene_symbol[TempIndices[!is.na(TempIndices)]]
UniqueGeneSymbols <- paste(NonUniqueGeneSymbols, 1:length(NonUniqueGeneSymbols), sep="_")
row.names(raw_data) <- UniqueGeneSymbols
colnames(raw_data) <- pheno_data$X

row.names(pheno_data) <- as.character(pheno_data$X)
pheno_data <- pheno_data[,-1]

# Finally, wrap this matrix up in a Seurat Object
data <- CreateSeuratObject(counts=raw_data,
                           project="basic_analysis",
                           min.cells=3,
                           min.features=200,
                           names.delim=NULL,
                           meta.data = pheno_data)

# First, find all mitochondrial genes, and count them as a percentage of total reads/cell
# In mouse, mitochondrial genes start with "mt-" so find all genes that match that pattern
# If you were doing this in a human dataset the pattern would be "^MT-"
data[["percent_mt"]] <- PercentageFeatureSet(object=data, pattern="^mt-")


# Typically, you would use much lower thresholds for mitochondrial genes (< 5%)
# This data set has lots of highly expressed mitochondrial genes though, so we'll leave them
quantnCountRNA <- quantile(data@meta.data$nCount_RNA, 0.05)
data <- subset(x=data, subset=nFeature_RNA > 200 & nCount_RNA > quantnCountRNA & percent_mt < 20)

print(sprintf("After filtering outliers: %d cells and %d genes", ncol(data), nrow(data)))

# # For raw count data, we would typically do LogNormalization:
# data <- NormalizeData(object=data, normalization.method="LogNormalize", scale.factor=10000)
# Again, these are the defaults, generate 2000 features using the "vst" feature selection method
# data <- FindVariableFeatures(object=data, selection.method="vst", nfeatures=2000)
# 
# 
# # Rescale all the genes
# scale_genes <- rownames(data)
# # If this takes too long, you can only rescale the variable genes
# # scale_genes <- VariableFeatures(object=data)
# data <- ScaleData(object=data, features=scale_genes)
# 
# # Use the highly variable genes to find principal components
# data <- RunPCA(object=data, features=VariableFeatures(object=data))
# 
# data <- RunTSNE(object=data, dims=1:15)
# data <- FindNeighbors(data, dims = 1:15)
# data <- FindClusters(object = data, resolution = 0.5)
# DimPlot(data, label = TRUE, reduction = "tsne")
# DimPlot(data, reduction = "tsne", group.by = "Individual")
# DimPlot(data, reduction = "tsne", group.by = "Time")

data <- SCTransform(data, method="qpoisson", vars.to.regress = NULL)
data <- RunPCA(data, verbose = FALSE)
data <- RunTSNE(data, dims = 1:30, verbose = FALSE)

data <- FindNeighbors(data, dims = 1:30, verbose = FALSE)
data <- FindClusters(data, verbose = FALSE)
DimPlot(data, label = TRUE, reduction = "tsne") 
DimPlot(data, reduction = "tsne", group.by = "Individual")
DimPlot(data, reduction = "tsne", group.by = "Time")
DimPlot(data, reduction = "tsne", group.by = "batch")
FeatureScatter(object=data, feature1="nFeature_RNA", feature2="PC_1")
FeaturePlot(object = data, features = "nFeature_RNA")
FeaturePlot(object = data, features = "PC_1")

# MarkersRes <- FindAllMarkers(data, assay = "SCT", slot = "data", test.use = "wilcox", return.thresh = 1e-6)
# MarkersRes1 <- FindMarkers(data, ident.1 = 0, assay = "SCT",slot = "data", test.use = "wilcox", return.thresh = 1e-6)
cluster0_data <- subset(x=data, subset=(seurat_clusters==0))
MarkersRes0 <- FindMarkers(cluster0_data, ident.1 = "5 day", group.by = "Time", assay = "SCT",slot = "data", test.use = "wilcox", return.thresh = 1e-6)
VlnPlot(cluster0_data, features = "Tac1-17705", group.by = "Time")
VlnPlot(cluster0_data, features = "Gm10116-5887", group.by = "Time")
VlnPlot(cluster0_data, features = "Hspa1l-9322", group.by = "Time")

# pb@assays@data[[1]][row.names(pb@assays@data[[1]]) == "Tac1-17705",]
# MarkersRes0[row.names(MarkersRes0)=="Tac1-17705",]
# pb@assays@data[[1]][row.names(pb@assays@data[[1]]) == "Gm10116-5887",]
# temp_ClusterResult <- tbl[[1]]
# temp_topClusterResult <- temp_ClusterResult[temp_ClusterResult$p_adj.loc < 1, ]
# temp_topClusterResult <- temp_topClusterResult[order(temp_topClusterResult$p_adj.loc),]
# temp_topClusterResult <- data.frame(Gene=row.names(temp_topClusterResult), temp_topClusterResult)
# temp_topClusterResult[row.names(temp_topClusterResult)=="Gm10116-5887", ]
# 
# pb@assays@data[[1]][row.names(pb@assays@data[[1]]) == "Hspa1l-9322",]
# temp_ClusterResult <- tbl[[1]]
# temp_topClusterResult <- temp_ClusterResult[temp_ClusterResult$p_adj.loc < 1, ]
# temp_topClusterResult <- temp_topClusterResult[order(temp_topClusterResult$p_adj.loc),]
# temp_topClusterResult <- data.frame(Gene=row.names(temp_topClusterResult), temp_topClusterResult)
# temp_topClusterResult[row.names(temp_topClusterResult)=="Hspa1l-9322", ]
# 
cluster1_data <- subset(x=data, subset=(seurat_clusters==1))
MarkersRes1 <- FindMarkers(cluster1_data, ident.1 = "5 day", group.by = "Time", assay = "SCT",slot = "data", test.use = "wilcox", return.thresh = 1e-6)
VlnPlot(cluster1_data, features = "Il1rl1-185", group.by = "Time")
# pb@assays@data[[2]][row.names(pb@assays@data[[1]]) == "Il1rl1-185",]
# MarkersRes1[row.names(MarkersRes1)=="Il1rl1-185",]


##differential analyses
data[["sTime"]] <- (data@meta.data$Time) %>%
  as.character() %>%
    gsub(" ", "_", .) %>%
    make.names()
data[["sIndividual"]] <- (data@meta.data$Individual) %>%
  as.character() %>%
  paste0("Individual_",.)

PhenoData <- data@meta.data

sce <- SummarizedExperiment(assays=list(counts=data@assays$RNA@counts, logcounts=data@assays$RNA@data), colData=PhenoData)
sce <- as(sce, "SingleCellExperiment")


(sce <- prepSCE(sce, 
                cluster_id = "seurat_clusters", # subpopulation assignments
                group_id = "sTime",  # group IDs
                sample_id = "sIndividual",   # sample IDs 
                drop = TRUE))  # drop all other colData columns

nk <- length(kids <- levels(sce$cluster_id))
ns <- length(sids <- levels(sce$sample_id))
names(kids) <- kids; names(sids) <- sids

toKeep <- table(sce$cluster_id, sce$sample_id) %>% 
  t() %>% 
  apply(., 2, function(x) median(x)) %>% 
  subset(., is_greater_than(., 5)) %>% 
  names()
sce <- subset(sce, , cluster_id %in% toKeep)

pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

(pb_mds <- pbMDS(pb))

##set up the design matrix
group_id <- colData(pb)[,1]
mm <- model.matrix(~ group_id)
colnames(mm) <- levels(pb$group_id)
# run DS analysis
res <- pbDS_update(pb, design=mm, coef=c(2),method="edgeR", verbose=TRUE, min_cells = 3)
# run DS analysis
# res <- pbDS(pb,method="edgeR", verbose=TRUE, min_cells = 3)

tbl <- res$table
ClusterNo <- 0
p_thresh <- 0.05
topClusterResults <- NULL
for(i in 1:length(tbl)) {
  temp_ClusterResult <- tbl[[i]]
  temp_topClusterResult <- temp_ClusterResult[temp_ClusterResult$p_adj.loc < p_thresh, ]
  temp_topClusterResult <- temp_topClusterResult[order(temp_topClusterResult$p_adj.loc),]
  temp_topClusterResult <- data.frame(Gene=row.names(temp_topClusterResult), temp_topClusterResult)
  topClusterResults <- rbind(topClusterResults, temp_topClusterResult)
}

##between-cluster comparsions
Ncells <- as.data.frame(metadata(pb)$n_cells)
Ncells <- Ncells %>%
  filter(Freq > 0)
short_Ncells <- Ncells %>% 
  spread(sample_id, Freq)

TotalCells <- colSums(short_Ncells[,-1], na.rm = T)
names(TotalCells) <- colnames(short_Ncells)[-1]
Clusters <- unique(Ncells$cluster_id)
SampleInfo <- colData(pb)
estimateCellStateChange <- function(k, Ncells, TotalCells, SampleInfo) {
  require(lme4)
  require(gdata)
  print(paste("Cluster", k))
  Ncells_sub <- Ncells %>%
    filter(cluster_id==k)
  FitData <- NULL
  TempIndices <- match(Ncells_sub$sample_id, names(TotalCells))
  TotalCells_sub <- TotalCells[TempIndices]
  for(i in 1:length(TotalCells_sub)) {
    InCluster=c(rep(1, Ncells_sub$Freq[i]), rep(0, (TotalCells_sub[i]-Ncells_sub$Freq[i])))
    Time=rep(SampleInfo$group_id, TotalCells_sub[i])
    ID=rep(row.names(SampleInfo)[i], TotalCells_sub[i])
    TempData <- data.frame(ID, InCluster, Time)
    FitData <- rbind(FitData, TempData)
  }
  FitData$Time <- relevel(FitData$Time, ref="X5_day")
  FitData$InCluster <- as.factor(FitData$InCluster)
  glmerFit1 <- glmer(InCluster ~ (1|ID) + Time, data=FitData, family = "binomial")
  sglmerFit1 <- summary(glmerFit1)
  TempRes1 <- (sglmerFit1$coefficients[-1,])


  return(TempRes1)
}

ClusterRes <- sapply(Clusters, estimateCellStateChange, Ncells, TotalCells, SampleInfo)
ClusterRes %<>% 
  as.data.frame() %>% 
  t() 
row.names(ClusterRes) <-  paste0("Cluster", Clusters)
ClusterRes <- data.frame(ClusterRes)
colnames(ClusterRes)[c(1,4)] <- c("logOddsRatio_11day_vs_5day","pvalue")
ClusterRes <- data.frame(ClusterRes, p.adjust= p.adjust(ClusterRes$pvalue, method = "BH"))

head(ClusterRes)

##batch correction
batch.list <- SplitObject(data, split.by = "batch")

for (i in 1:length(batch.list)) {
  batch.list[[i]] <- SCTransform(batch.list[[i]], verbose = FALSE, method="qpoisson", vars.to.regress = NULL)
}

batch.features <- SelectIntegrationFeatures(object.list = batch.list, nfeatures = 3000)
batch.list <- PrepSCTIntegration(object.list = batch.list, anchor.features = batch.features, 
                                    verbose = FALSE)
batch.anchors <- FindIntegrationAnchors(object.list = batch.list, normalization.method = "SCT", 
                                           anchor.features = batch.features, verbose = FALSE, k.filter  = 10)
batch.integrated <- IntegrateData(anchorset = batch.anchors, normalization.method = "SCT", 
                                     verbose = FALSE)

batch.integrated <- RunPCA(batch.integrated, verbose = FALSE)
batch.integrated <- RunTSNE(batch.integrated, dims = 1:30, verbose = FALSE)

batch.integrated <- FindNeighbors(batch.integrated, dims = 1:30, verbose = FALSE)
batch.integrated <- FindClusters(batch.integrated, verbose = FALSE)
DimPlot(batch.integrated, label = TRUE, reduction = "tsne") 
DimPlot(batch.integrated, reduction = "tsne", group.by = "Individual")
DimPlot(batch.integrated, reduction = "tsne", group.by = "Time")
DimPlot(batch.integrated, reduction = "tsne", group.by = "batch")

