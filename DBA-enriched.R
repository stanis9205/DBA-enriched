setwd("/local/storage/JQ1/JQ1-R-Outputs/3W-vs-2M/")

# Loading Packages ---------------------------------------------------------

set.seed(1234)

library(Seurat)
library(dplyr)
library(Matrix)
library(sva)
library(SingleR)
library(reshape2)
library(pheatmap)
library(kableExtra)
library(biomaRt)
library(celda)
library(scDblFinder)
library(SoupX) # install.packages("SoupX")
library(tidyverse) # install.packages("tidyverse")
library(knitr)

if(!dir.exists("output")) dir.create("output")
if(!dir.exists("data")) dir.create("data")

# install.packages("DropletUtils")


# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("scDblFinder")

# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("celda")



########################################################################
#
#  1 Data preprocessing
# 
# ######################################################################
#======Load the data files and Set up Seurat object =========================
# Load the dataset

# setup Seurat objects since both count matrices have already filtered
# cells, we do no additional filtering here

df_samples <- readxl::read_excel("sample_info.xlsx")
(samples <- df_samples$Samples)
(conditions <- df_samples$Conditions)
(experiments <- df_samples$Experiments)

# testis_raw <- list()
# testis_Seurat <- list()


for(i in 1:length(samples)){
  testis_raw[[i]] <- Read10X(data.dir = paste0("/local/storage/JQ1/JQ1-R-Outputs/3W-vs-2M/data/",
                                               samples[i],"/outs/filtered_feature_bc_matrix/"))
  colnames(testis_raw[[i]]) <- paste0(samples[i],"_",colnames(testis_raw[[i]]))
  testis_Seurat[[i]] <- CreateSeuratObject(testis_raw[[i]],
                                           min.cells = 3,
                                           min.genes = 200,
                                           names.delim = "_",
                                           project = "3W-vs-2M")
  # These next lines of code add different information (lines 28 & 29) to the meta-data
  testis_Seurat[[i]]@meta.data$conditions <- conditions[i]
  testis_Seurat[[i]]@meta.data$experiments <- experiments[i]
}


# testis <- Reduce(function(x, y) merge(x, y), testis_Seurat)

# remove(testis_raw,testis_Seurat);gc()

# testis <- subset(testis, subset = nFeature_RNA > 200)

# save(testis, file = "./data/testis_loaded.Rda")
# load(file = "./data/testis_loaded.Rda")

################################################################################
#======1.2 QC, pre-processing and normalizing the data=========================
################################################################################

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
genes_all <- getBM(attributes=c("chromosome_name", "start_position","end_position","external_gene_name"), 
                   filters="external_gene_name",
                   values=rownames(testis@assays$RNA@counts),
                   mart = ensembl)

genes_y <- genes_all[genes_all$chromosome_name=="Y",]
genes_X <- genes_all[genes_all$chromosome_name=="X",]

testis[["percent.mt"]] <- PercentageFeatureSet(testis, pattern = "^mt-")
testis[["percent.XY"]] <- PercentageFeatureSet(testis, features = genes_all$external_gene_name[genes_all$chromosome_name %in% c("X", "Y")])

VlnPlot(testis, features = c("nFeature_RNA"), ncol = 1, pt.size=0)
VlnPlot(testis, features = c("nCount_RNA"), ncol = 1, pt.size=0)
VlnPlot(testis, features = c("percent.mt"), ncol = 1, pt.size=0)

plot1 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


testis <- subset(testis, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & nCount_RNA > 1000 & percent.mt < 5)


VlnPlot(testis, features = c("nFeature_RNA"), ncol = 1, group.by = "orig.ident", pt.size=0)
VlnPlot(testis, features = c("nCount_RNA"), ncol = 1, group.by = "orig.ident", pt.size=0)
VlnPlot(testis, features = c("percent.mt"), ncol = 1, group.by = "orig.ident", pt.size=0)

plot1 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(testis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

################################################################################
#==============scDblFinder======================================================
################################################################################

# Cell doublets will be identified using scDblFinder 
# (https://bioconductor.org/packages/release/bioc/vignettes/scDblFinder/inst/doc/introduction.html)
# which uses functions from various other doublet detection methods, but has been
# rewritten to be more compute efficient. A comparison between other methods
# and scDblFinder can be found in [Xi and Li](https://arxiv.org/pdf/2101.08860.pdf)
# In this method, doublets are detected through the introduction of pseudo-doublets
# by combining two random cells. A nearest-neighbor approach is then employed to 
# see which cells associate most closely with these pseudo-doublest. The percent
# of cells called as doublets is based on 10x genomics known doublet rate.

set.seed(1234)

# testis[["samples"]] <- testis$orig.ident

filt.matrix <-list()
raw.matrix <-list()

for(i in 1:length(samples)){
  filt.matrix[[i]] <- Read10X(data.dir = paste0("/local/storage/JQ1/JQ1-R-Outputs/3W-vs-2M/data/",
                                                samples[i],"/outs/filtered_feature_bc_matrix/"))
  raw.matrix[[i]] <- Read10X(data.dir = paste0("/local/storage/JQ1/JQ1-R-Outputs/3W-vs-2M/data/",
                                               samples[i],"/outs/raw_feature_bc_matrix/"))
  colnames(filt.matrix[[i]]) <- paste0(samples[i],"_",colnames(filt.matrix[[i]]))
  colnames(raw.matrix[[i]]) <- paste0(samples[i],"_",colnames(raw.matrix[[i]]))
}







sce <- CreateSeuratObject(counts = filt.matrix)

# sce$colData <- filt.matrix[["samples"]]
sce <- scDblFinder(sce, samples="samples")

sce_sub <- subset(sce, , scDblFinder.class=="singlet")
singlet_cell_id <- colnames(sce_sub)

table(sce$scDblFinder.class)
table(sce_sub$scDblFinder.class)


sc <- SoupChannel(raw.matrix, filt.matrix)

#Convert to single cell experiment
sce <- as.SingleCellExperiment(filt.matrix, raw.matrix)

# Run scDblFinder
sce$colData <- testis[["samples"]]
sce <- scDblFinder(sce, samples="samples")

sce_sub <- subset(sce, , scDblFinder.class=="singlet")
singlet_cell_id <- colnames(sce_sub)

table(sce$scDblFinder.class)
table(sce_sub$scDblFinder.class)

# For this dataset, it appears that 176414 cells are singlets and 15142 are doublets

# I will then use the sce_sub I have created, as I want to keep only the singlet
# cells for further downstream use.
testis <- as.Seurat(sce_sub, counts = "counts", data = "logcounts")

save(sce_sub, file = "./data/sce_sub.Rda")

################################################################################
#======================SoupX==================================================
################################################################################

# From the Sethupathy lab: In single cell experiments, non-endogenous RNA can 
# make up anywhere from 2-50% of all counts, although 10-20% seems like a decent 
# expectation. To remove these exogenous counts, we can use the SoupX program that 
# models the background expression from empty droplets. More information can be 
# found in the [SoupX manuscript]
# (https://academic.oup.com/gigascience/article/9/12/giaa151/6049831) or the 
# [SoupX tutorial](https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html). 
# I further took the advice of this [CellGen vignette]
# (https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html) to use the 
# Seurat clustering. Other options for ambient RNA removal include Cell Bender, 
# DecontX, and SoupOrCell.

# According to my conversations with Adriana, the Yao lab thinks SoupX might be 
# overcorrecting. Therefore, they advised me to use SoupX to get a contamination
# estimate, then use those estimates in the DeContX pipeline. 

sce_sub_SoupX <- runSoupX(sce_sub, sample="samples")






# I am starting here with getting the basics of DeContX
sce_sub <- decontX(x = sce_sub, batch= "samples")

plotDecontXContamination(sce_sub)

# After removing unwanted cells from the dataset, the next step is to integrate
# The first way is following the Seurat vignette. The second way is to split the dataset
# into different Seurat objects based on 3W or 2M or UC to integrate biological 
# replicates. Then will merge.

################################################################################
#=============Large Dataset Integration Vignette================================
################################################################################

# This will integrate everyone together, using the three UCs as references
# https://satijalab.org/seurat/articles/integration_large_datasets.html

# In this vignette, we present a slightly modified workflow for the integration 
# of scRNA-seq datasets. Instead of utilizing canonical correlation analysis (‘CCA’) 
# to identify anchors, we instead utilize reciprocal PCA (‘RPCA’). When determining 
# anchors between any two datasets using RPCA, we project each dataset into the 
# others PCA space and constrain the anchors by the same mutual neighborhood 
# requirement. The commands for both workflows are largely identical, but the two 
# methods may be applied in different context.

# By identifying shared sources of variation between datasets, CCA is well-suited 
# for identifying anchors when cell types are conserved, but there are very 
# substantial differences in gene expression across experiments. CCA-based 
# integration therefore enables integrative analysis when experimental conditions 
# or disease states introduce very strong expression shifts, or when integrating 
# datasets across modalities and species. However, CCA-based integration may also 
# lead to overcorrection, especially when a large proportion of cells are 
# non-overlapping across datasets.

# RPCA-based integration runs significantly faster, and also represents a more 
# conservative approach where cells in different biological states are less 
# likely to ‘align’ after integration. We therefore,recommend RPCA during 
# integrative analysis where: * A substantial fraction of cells in one dataset 
# have no matching type in the other * Datasets originate from the same platform 
# (i.e. multiple lanes of 10x genomics) * There are a large number of datasets 
# or cells to integrate (see INSERT LINK for more tips on integrating large datasets)

testis.list.V <- SplitObject(testis, split.by = "orig.ident")

# normalize and identify variable features for each dataset independently

testis.list.V <- lapply(X = testis.list.V, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across data sets for integration
# run PCA on each dataset using these features
features <- SelectIntegrationFeatures(object.list = testis.list.V)

testis.list.V <- lapply(X = testis.list.V, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

JQ1.anchors.V <- FindIntegrationAnchors(object.list = testis.list.V, reference =c(25,26,27), 
                                        reduction = "rpca", anchor.features = features)

# this command creates an 'integrated' data assay
JQ1.combined.V <- IntegrateData(anchorset = JQ1.anchors.V, dims = 1:50)

# specify that we will perform downstream analysis on the corrected data 
# note that the original unmodified data still resides in the 'RNA' assay
DefaultAssay(JQ1.combined.V) <- "integrated"

# Run the standard workflow for visualization and clustering
JQ1.combined.V <- ScaleData(JQ1.combined.V, verbose = FALSE)
JQ1.combined.V <- RunPCA(JQ1.combined.V, npcs = 30, verbose = FALSE)
JQ1.combined.V <- RunUMAP(JQ1.combined.V, reduction = "pca", dims = 1:30)

# Reslolution 3
JQ1.combined.V <- FindNeighbors(JQ1.combined.V, reduction = "pca", force.recalc = T, k.param = 30, dims = 1:30) %>% 
  FindClusters(resolution = 3, n.start = 100)
JQ1.combined.V[["res.3"]] <- Idents(object = JQ1.combined.V)

# Resolution 2
JQ1.combined.V <- FindNeighbors(JQ1.combined.V, reduction = "pca", force.recalc = T, k.param = 30, dims = 1:30) %>% 
  FindClusters(resolution = 2, n.start = 100)
JQ1.combined.V[["res.2"]] <- Idents(object = JQ1.combined.V)

# Resolution 1.5
JQ1.combined.V <- FindNeighbors(JQ1.combined.V, reduction = "pca", force.recalc = T, k.param = 30, dims = 1:30) %>% 
  FindClusters(resolution = 1.5, n.start = 100)
JQ1.combined.V[["res.1.5"]] <- Idents(object = JQ1.combined.V)

# Visualization
p1 <- DimPlot(JQ1.combined.V, reduction = "umap", group.by = "conditions")
p2 <- DimPlot(JQ1.combined.V, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + NoAxes()
# p2 + NoAxes()
p1 + p2

p3 <- DimPlot(JQ1.combined.V, reduction = "umap", split.by = "conditions")
p3 + NoAxes()

# 3W Only
JQ1.combined.V.3W <- JQ1.combined.V[, JQ1.combined.V@meta.data[, "experiments"] == "3W"]
p1 <- DimPlot(JQ1.combined.V.3W, reduction = "umap", group.by = "conditions", split.by = "conditions")
p2 <- DimPlot(JQ1.combined.V.3W, reduction = "umap", label = TRUE, repel = TRUE)
p1 + NoAxes()
p2 <- p2 + NoAxes()
p1

# 2M Only
JQ1.combined.V.2M <- JQ1.combined.V[, JQ1.combined.V@meta.data[, "experiments"] == "2M"]
p1 <- DimPlot(JQ1.combined.V.2M, reduction = "umap", split.by = "conditions")
p2 <- DimPlot(JQ1.combined.V.2M, reduction = "umap", label = TRUE, repel = TRUE)
p1 <- p1 + NoAxes()
p2 <- p2 + NoAxes()
p1

# UC Only
JQ1.combined.V.UC <- JQ1.combined.V[, JQ1.combined.V@meta.data[, "experiments"] == "UC"]
p1 <- DimPlot(JQ1.combined.V.UC, reduction = "umap", group.by = "conditions")
p2 <- DimPlot(JQ1.combined.V.UC, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2

save(JQ1.combined.V, file = "./data/JQ1.combined.V.Rda")

# load("./data/JQ1.combined.V.Rda")


#============== Identifying cell types =======================

FeaturePlot(JQ1.combined.V, feature = "Ptgds", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

# Hematopoietic & Endothelial cell marker
FeaturePlot(JQ1.combined.V, feature = "Laptm5", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Hbb-bt", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Vcam1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Insl3", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Hsd3b1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Fabp3", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

HandE_mark <- c("Laptm5", "Hbb-bt", "Vcam1","Insl3","Hsd3b1","Fabp3")
FeaturePlot(JQ1.combined.V, features = HandE_mark, label = F, reduction = "umap", min.cutoff = 0)

# Smooth muscle marker
FeaturePlot(JQ1.combined.V, feature = "Acta2", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Col1a2", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

smooth_mark <- c("Acta2", "Col1a2")
FeaturePlot(JQ1.combined.V, features = smooth_mark, label = F, reduction = "umap", min.cutoff = 0)

# Sertoli cell marker
FeaturePlot(JQ1.combined.V, feature = "Cst9", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Cldn11", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Mro", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

sertoli_mark <- c("Cst9", "Cldn11", "Mro")
FeaturePlot(JQ1.combined.V, features = sertoli_mark, label = F, reduction = "umap", min.cutoff = 0)

# Undifferentiated SG
FeaturePlot(JQ1.combined.V, feature = "Sall4", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Zbtb16", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Gfra1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Dazl", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Kit", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Neurog3", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()


# Differentiated SG
FeaturePlot(JQ1.combined.V, feature = "Sohlh1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Stra8", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Gm4969", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Prdm9", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()


# Early Spermatocytes
FeaturePlot(JQ1.combined.V, feature = "Fmr1nb", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Prss50", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Defb19", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()


# Scytes
FeaturePlot(JQ1.combined.V, feature = "Hormad1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Sycp1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

# Round spermatid
FeaturePlot(JQ1.combined.V, feature = "Lyzl1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Acrv1", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()
FeaturePlot(JQ1.combined.V, feature = "Hemgn", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

# Spermatid
FeaturePlot(JQ1.combined.V, feature = "Tssk6", label = F, reduction = "umap", min.cutoff = 0) + NoAxes()

germ_mark <- c("Gfra1","Zbtb16","Sall4","Dmrt1","Dazl","Kit","Cdca8",
               "Id4","Sycp3","Tesmin","Nxt1","Shcbp1l","Aurka","Lyzl1",
               "Acrv1","Hemgn","Txndc8","Tssk6","Oaz3","Prm2")

FeaturePlot(JQ1.combined, feature = germ_mark[c(3:5,8:9,12,15,19:20)], label = F, reduction = "umap", raster=FALSE)


# kable(t(table(testis@meta.data[,c(1,19)]))) %>% kableExtra::kable_styling()

JQ1.combined.V[["res.3_ordered"]] <- factor(unlist(JQ1.combined.V[["res.3"]]), levels = rev(c(14,32,26,9,16,17,19,10,12,22,18,8,11,6,0,7,3,4,15,33,23,25,24,13,5,29,2,21,30)), ordered = T)

DotPlot(JQ1.combined.V, features = germ_mark,
        cols = c("blue","red"), dot.scale = 8, group.by = "res.3_ordered")

########################################################

# Resolution 1.5 Cell Annotation

JQ1.combined.V[["res.1.5_cells"]] <- as.character(as.numeric(JQ1.combined.V@meta.data$res.1.5)-1)

JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 20] <- "Smooth_Muscle"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 38] <- "Smooth_Muscle"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 31] <- "Hematopoetic_Cells"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 39] <- "Hematopoetic_Cells"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 33] <- "Red_Blood_Cells"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 30] <- "Endothelial_Cells"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 29] <- "Sertoli_Cells"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 5] <- "Sertoli_Cells"

JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 36] <- "Undifferentiated_Spermatogonia"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 35] <- "Differentiated_Spermatogonia"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 15] <- "Differentiated_Spermatogonia"

JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 27] <- "Early_Spermatocytes"

JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 17] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 18] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 21] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 34] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 25] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 11] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 24] <- "Spermatocytes"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 23] <- "Spermatocytes"

JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 2] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 8] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 6] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 13] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 9] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 3] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 14] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 37] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 22] <- "Spermatids"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 28] <- "Spermatids"

JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 16] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 7] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 19] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 26] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 32] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 12] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 1] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 0] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 10] <- "Spermatozoa"
JQ1.combined.V[["res.1.5_cells"]][JQ1.combined.V[["res.1.5_cells"]] == 4] <- "Spermatozoa"

########################################################

# Resolution 2 Cell Annotation

JQ1.combined.V[["res.2_cells"]] <- as.character(as.numeric(JQ1.combined.V@meta.data$res.2)-1)

JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 18] <- "Smooth_Muscle"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 49] <- "Smooth_Muscle"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 40] <- "Hematopoetic_Cells"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 46] <- "Hematopoetic_Cells"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 44] <- "Red_Blood_Cells"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 36] <- "Endothelial_Cells"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 30] <- "Sertoli_Cells"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 2] <- "Sertoli_Cells"

JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 47] <- "Undifferentiated_Spermatogonia"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 12] <- "Differentiated_Spermatogonia"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 45] <- "Differentiated_Spermatogonia"

JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 35] <- "Early_Spermatocytes"

JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 21] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 26] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 25] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 37] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 23] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 38] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 39] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 43] <- "Spermatocytes"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 5] <- "Spermatocytes"

JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 8] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 41] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 16] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 7] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 4] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 10] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 6] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 19] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 31] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 27] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 15] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 34] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 32] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 22] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 48] <- "Spermatids"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 29] <- "Spermatids"


JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 28] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 33] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 24] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 14] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 17] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 9] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 13] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 42] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 11] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 20] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 1] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 0] <- "Spermatozoa"
JQ1.combined.V[["res.2_cells"]][JQ1.combined.V[["res.2_cells"]] == 3] <- "Spermatozoa"


save(JQ1.combined.V, file = "output/JQ1.combined.V_clustered_with_identities.Rda")








testis[["res.4_cells_genotype"]] <- testis[["res.4_cells"]]

testis[["res.4_cells_genotype"]][testis[["res.4_cells"]] == "Pre-Leptotene" & testis@meta.data$orig.ident == "WTh"] <- "WTh_Pre-Leptotene"
testis[["res.4_cells_genotype"]][testis[["res.4_cells"]] == "Pre-Leptotene" & testis@meta.data$orig.ident == "KEc"] <- "KEc_Pre-Leptotene"

testis[["res.4_cells_genotype"]][testis[["res.4_cells"]] == "Spermatocyte1" & testis@meta.data$orig.ident == "WTh"] <- "WTh_Spermatocyte1"
testis[["res.4_cells_genotype"]][testis[["res.4_cells"]] == "Spermatocyte2" & testis@meta.data$orig.ident == "WTh"] <- "WTh_Spermatocyte2"

testis[["res.4_cells_genotype"]][testis[["res.4_cells"]] == "Spermatocyte1" & testis@meta.data$orig.ident == "KEc"] <- "KEc_Spermatocyte1"
testis[["res.4_cells_genotype"]][testis[["res.4_cells"]] == "Spermatocyte2" & testis@meta.data$orig.ident == "KEc"] <- "KEc_Spermatocyte2"

unique(testis[["res.4_cells_genotype"]])

tmp <- testis@meta.data$res.3_cells_genotype
names(tmp) <- names(Idents(testis))
Idents(testis) <- tmp; rm(tmp);GC()

DimPlot(testis, reduction = "umap", label = T, label.size = 4)
DimPlot(testis, reduction = "umap", label = F, label.size = 4, group.by = "cell_idents_4")

testis[[i]]@meta.data$conditions <- experiments[i]

testis.subset <- subset(x = testis, subset = (conditions == c("Untreated-Control", "Treatment")))

DimPlot(object = testis, reduction = "umap", label = F, group.by = "cell_idents_4", cells.highlight = testis@meta.data$conditions == "Untreated-Control", cols.highlight = "darkred", cols= "grey")

DimPlot(testis, reduction = "umap", label = F, label.size = 4, raster=FALSE, group.by = "cell_idents_4")

DimPlot(testis.subset, reduction = "umap", label = F, label.size = 3, raster=FALSE, group.by = "conditions")

save(testis, file = "output/testis_clustered_with_identities.Rda")


################################################################################

KEc_spermatocytes2.markers_all_negbinom <- FindMarkers(testis, ident.1 = "KEc_Spermatocyte2", ident.2 = "WTh_Spermatocyte2", min.pct = 0, logfc.threshold = 0, test.use = "negbinom")
head(KEc_spermatocytes2.markers_all_negbinom)

KEc_spermatocytes2.markers_all_negbinom[rownames(KEc_spermatocytes2.markers_all_negbinom) %in% genes_y$external_gene_name,]

KEc_spermatocytes2.markers_all_negbinom[rownames(KEc_spermatocytes2.markers_all_negbinom) %in% genes_all$external_gene_name[genes_all$chromosome_name=="X"],]


