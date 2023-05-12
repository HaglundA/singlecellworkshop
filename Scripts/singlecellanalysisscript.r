


library(Seurat)
library(dplyr)

pbmc.data <- Read10X(data.dir = "/rds/general/user/ah3918/projects/bioinformatics_module/live/singlecell_tutorial/cellranger_output/filtered_gene_bc_matrices/hg19")
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]


pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)


### Create ViolinPlots and Scattrplots
plot1<-VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


pdf("Violinplot.pdf")
plot1 
dev.off()

pdf("Scatterplot.pdf")
plot2 
dev.off()


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
plot3<-VlnPlot(pbmc,features="percent.mt", ncol = 1)
pdf("Violinplot_Mitochondrialcontent.pdf")
plot3
dev.off()


pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pdf("Violinplot_postQC.pdf")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA","percent.mt"), ncol = 3)
dev.off()



######### Data normalization


message("Normalizing and scaling the data...")
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- NormalizeData(pbmc)

pdf("Histogram_normalized_data.pdf")
hist(as.matrix(pbmc[["RNA"]]@data))
hist(as.matrix(pbmc[["RNA"]]@data),xlim = c(0.5,4),ylim=c(0,1e6))
dev.off()

pdf("Histogram_nonnormalized_data.pdf")
hist(as.matrix(pbmc[["RNA"]]@counts))
hist(as.matrix(pbmc[["RNA"]]@counts),xlim = c(0.5,50),ylim=c(0,1e6),breaks=100)
dev.off()


message("Finding variable features..")

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


pdf("Variable_features.pdf")
plot2 
dev.off()


pbmc <- ScaleData(pbmc)


message("running PCA analysis")

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

pdf("PC_loadings.pdf")
VizDimLoadings(pbmc, dims = 1:4, reduction = "pca")
dev.off()


pdf("PC1_heatmap.pdf")
DimHeatmap(pbmc, dims = 1,cells=500,balanced=TRUE)
dev.off()

message("running jackstraw..")
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

pdf("JackstrawPlot.pdf")
JackStrawPlot(pbmc, dims = 1:15)
dev.off()

pdf("Elbowplot.pdf")
ElbowPlot(pbmc)
dev.off()




##### Clustering and annotation
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:10)
pbmc<-RunTSNE(pbmc,dims=1:10)


pdf("UMAP_vs_tsne.pdf")
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")
dev.off()


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


markers=pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)
markers=markers$gene

pdf("topmarkers_violinplot.pdf")
VlnPlot(pbmc, features = markers[1:2], slot = "counts", log = TRUE)
dev.off()

pdf("topmarkers_featureplot.pdf")
FeaturePlot(pbmc, features = markers[1:5])
dev.off()

pdf("markergenes_featureplot.pdf")
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
    "CD8A"))
dev.off()


pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10

pdf("topgenes_heatmap.pdf")
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()

pdf("MS4A1_example.pdf")
p1=DimPlot(pbmc, reduction = "umap")
p2=FeaturePlot(pbmc,features="MS4A1")
p1+p2
dev.off()


markergenes=c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP","CD8A")

pdf("markers_heatmap.pdf")
DoHeatmap(pbmc, features = markergenes)+NoLegend()
dev.off()

message("Annotating cells..")

levels(pbmc)
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids)=levels(pbmc)

pbmc <- RenameIdents(pbmc, new.cluster.ids)

dev.off("UMAP_annotated_celltypes.pdf")
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()