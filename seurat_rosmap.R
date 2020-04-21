library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(Matrix)

### https://community.rstudio.com/t/r-studio-crashes-when-loading-seurat-package/57931/6

#### we have to put genes.tsv and barcodes.tsv as similar format as 
#### ./filtered_gene_bc_matrices/hg19/[genes.tsv/barcodes.tsv] 

### load raw rosmap single cell data and create a seurat object
not_fil_rosmap_scRNA_data = Read10X(data.dir = "../ROSMAP_scRNA_data/Processed/hg38/")

# Initialize the Seurat object with the raw data
not_fil_rosmap_scRNA = CreateSeuratObject(counts = not_fil_rosmap_scRNA_data, 
                                          project = "rosmap_scRNA", 
                                          min.cells = 10, min.features = 500)

not_fil_rosmap_scRNA

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
not_fil_rosmap_scRNA[["percent.mt"]] = PercentageFeatureSet(not_fil_rosmap_scRNA, pattern = "^MT-")

# Visualize QC metrics as a violin plot
pdf("../ROSMAP_scRNA_data/Processed/plots/feature_count_rosmap.pdf")
VlnPlot(not_fil_rosmap_scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3, pt.size = 0)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("../ROSMAP_scRNA_data/Processed/plots/feature_feature_rel_rosmap.pdf",width = 12, height = 9)
plot1 <- FeatureScatter(not_fil_rosmap_scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(not_fil_rosmap_scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1, plot2))
plot1 + plot2
dev.off()

### decide boundary based on the plot feature_feature_rel_rosmap.pdf

not_fil_rosmap_scRNA = subset(not_fil_rosmap_scRNA, 
                               subset = nFeature_RNA > 500 & nFeature_RNA < 8500 & 
                                percent.mt > -Inf & percent.mt < 15)
not_fil_rosmap_scRNA = NormalizeData(not_fil_rosmap_scRNA, 
                                      normalization.method = "LogNormalize", 
                                      scale.factor = 10000)
not_fil_rosmap_scRNA = NormalizeData(not_fil_rosmap_scRNA)

not_fil_rosmap_scRNA = FindVariableFeatures(not_fil_rosmap_scRNA, 
                                             selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(not_fil_rosmap_scRNA), 10)

# plot variable features with and without labels
pdf("../ROSMAP_scRNA_data/Processed/plots/top_10_genes.pdf", width = 12, height=8)
plot1 <- VariableFeaturePlot(not_fil_rosmap_scRNA)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#CombinePlots(plots = list(plot1, plot2))
plot1 + plot2
dev.off()

all.genes = rownames(not_fil_rosmap_scRNA)
not_fil_rosmap_scRNA = ScaleData(not_fil_rosmap_scRNA, features = all.genes)
not_fil_rosmap_scRNA = RunPCA(not_fil_rosmap_scRNA, 
                              features = VariableFeatures(object = not_fil_rosmap_scRNA))
print(not_fil_rosmap_scRNA[["pca"]], dims = 1:5, nfeatures = 5)

pdf("../ROSMAP_scRNA_data/Processed/plots/pca_loading.pdf")
VizDimLoadings(not_fil_rosmap_scRNA, dims = 1:2, reduction = "pca")
dev.off()

pdf("../ROSMAP_scRNA_data/Processed/plots/pca_2D.pdf")
DimPlot(not_fil_rosmap_scRNA, reduction = "pca")
dev.off()

DimHeatmap(not_fil_rosmap_scRNA, dims = 1, cells = 500, balanced = TRUE)
pdf("../ROSMAP_scRNA_data/Processed/plots/heatmap_PC.pdf", width = 12, height = 7)
DimHeatmap(not_fil_rosmap_scRNA, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
not_fil_rosmap_scRNA = JackStraw(not_fil_rosmap_scRNA, num.replicate = 100)
not_fil_rosmap_scRNA = ScoreJackStraw(not_fil_rosmap_scRNA, dims = 1:20)
JackStrawPlot(not_fil_rosmap_scRNA, dims = 1:15)
ElbowPlot(not_fil_rosmap_scRNA)

not_fil_rosmap_scRNA = FindNeighbors(not_fil_rosmap_scRNA, dims = 1:10)
not_fil_rosmap_scRNA = FindClusters(not_fil_rosmap_scRNA, resolution = 0.5)
head(Idents(not_fil_rosmap_scRNA), 5)

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
not_fil_rosmap_scRNA = RunUMAP(not_fil_rosmap_scRNA, dims = 1:10)

nbt=run_tsne(not_fil_rosmap_scRNA,dims.use = 1:10,max_iter=2000)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pdf("../ROSMAP_scRNA_data/Processed/plots/umap_rosmap.pdf")
DimPlot(not_fil_rosmap_scRNA, reduction = "umap")
dev.off()


# find all markers of cluster 1
cluster1.markers = FindMarkers(not_fil_rosmap_scRNA, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(not_fil_rosmap_scRNA, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells, report only the positive ones
not_fil_rosmap_scRNA.markers = FindAllMarkers(not_fil_rosmap_scRNA, 
                                               only.pos = TRUE, 
                                               min.pct = 0.25, logfc.threshold = 0.25)
not_fil_rosmap_scRNA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)


cluster1.markers = FindMarkers(not_fil_rosmap_scRNA, ident.1 = 0, 
                                logfc.threshold = 0.25, test.use = "roc", 
                                only.pos = TRUE)


#############
VlnPlot(not_fil_rosmap_scRNA, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(not_fil_rosmap_scRNA, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(not_fil_rosmap_scRNA, 
            features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
top10 = not_fil_rosmap_scRNA.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(not_fil_rosmap_scRNA, features = top10$gene) + NoLegend()

new.cluster.ids = c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) = levels(not_fil_rosmap_scRNA)
not_fil_rosmap_scRNA = RenameIdents(not_fil_rosmap_scRNA, new.cluster.ids)
DimPlot(not_fil_rosmap_scRNA, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


