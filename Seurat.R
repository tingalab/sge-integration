library(Seurat)
library(dplyr)
library(tidyverse)
library(BayesSpace)
library(purrr)
library(janitor)
library(magrittr)
library(purrr)
library(patchwork)

## ----- Load Normal samples ------
s1 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/BL_190222U/outs/") %>% 
  PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
  # NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  # FindVariableFeatures(nfeatures = 3000) %>% 
  # ScaleData(vars.to.regress = "percent.mt") %>%
  SCTransform(assay="Spatial", vars.to.regress = "percent.mt", return.only.var.genes = FALSE) %>%
  RunPCA()
s1$orig.ident <- "BL_190222U"

s2 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/BL_191561U/outs/") %>%
  PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
  # NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  # FindVariableFeatures(nfeatures = 3000) %>% 
  # ScaleData(vars.to.regress = "percent.mt") %>%
  SCTransform(assay="Spatial", vars.to.regress = "percent.mt", return.only.var.genes = FALSE) %>%
  RunPCA()
s2$orig.ident <- "BL_191561U"

s3 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/BL_197129U/outs/") %>% 
  PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
  # NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  # FindVariableFeatures(nfeatures = 3000) %>%
  # ScaleData(vars.to.regress = "percent.mt", features = row.names(s3)) %>%
  SCTransform(assay="Spatial", vars.to.regress = "percent.mt", return.only.var.genes = FALSE) %>%
  RunPCA()
s3$orig.ident <- "BL_197129U"

s4 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/Donor4_ureter/outs/") %>% 
  PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
  # NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  # FindVariableFeatures(nfeatures = 3000) %>% 
  # ScaleData(vars.to.regress = "percent.mt") %>%
  SCTransform(assay="Spatial", vars.to.regress = "percent.mt", return.only.var.genes = FALSE) %>%
  RunPCA()
s4$orig.ident <- "Donor4_ureter"

UMAP <- function(sample, n_dims=15, res=0.8){
  sample <- FindNeighbors(sample,reduction="pca", dims=1:n_dims)
  sample <- FindClusters(sample, resolution=res)
  sample <- RunUMAP(sample,reduction="pca", dims=1:n_dims)
  return(sample)
}

s.obj.list <- list(s1,s3)
features <- SelectIntegrationFeatures(s.obj.list , nfeatures=3000)
s.obj.list <- PrepSCTIntegration(s.obj.list, anchor.features=features)
# all.genes <- unique(c(row.names(1), row.names(s2), row.names(s3), row.names(s4)))
# all.genes <- all.genes[-grep(x = all.genes, pattern= "^MT-")]
sample.anchors <- FindIntegrationAnchors(s.obj.list, dims=1:30, reduction="cca", normalization.method="SCT", anchor.features=features)
#sample.anchors <- FindIntegrationAnchors(s.obj.list, dims=1:30, reduction="cca", scale = FALSE, normalization.method="LogNormalize", anchor.features=features)
s.obj.integrated <- IntegrateData(anchorset=sample.anchors, dims=1:30, normalization.method="SCT") %>%
  # ScaleData() %>%
  RunPCA() %>%
  UMAP(n_dims = 20, res=0.3)

SpatialDimPlot(s.obj.integrated, label=FALSE) %T>% ggsave(filename="ureter-integratedVisium.pdf", plot=., width=28,height=7,units='in') %>% plot
saveRDS(s.obj.integrated, file="integrated-ureter-visium.Rds")
tab<-table(Idents(s.obj.integrated),s.obj.integrated@meta.data$orig.ident)

#Convert to data frame
tab<-as.data.frame.matrix(tab)
tab<-cbind(rownames(tab),tab)
colnames(tab)[1]<-'cluster'

#Add totals

tab2<-adorn_totals(tab, c("row","col"))
write.csv (tab2,"ureter-cells_by_cluster_by_sample.csv",row.names = FALSE)
rm(tab)
rm(tab2)

# Find markers for found clusters
top_markers <- FindAllMarkers(s.obj.integrated, assay='integrated', group.by= "subclass") %>% 
  filter(p_val_adj<0.05) %>% 
  mutate(pct.diff=pct.1-pct.2) %>% 
  group_by(cluster) %>% 
  arrange(desc(pct.diff,avg_log2FC), cluster)

top10_pctdiff_markers <- top_markers %>%
  top_n(10,pct.diff) %>% 
  arrange(cluster)

top10_log2FC_markers <- top_markers %>%
  top_n(10,avg_log2FC) %>% 
  arrange(cluster)

DoHeatmap(s.obj.integrated, features=top10_pctdiff_markers$gene, slot='data', assay='integrated')
ggsave(filename="ureter-tumor-pctdiff-heatmap.pdf", width=16,height=10,units='in')
DoHeatmap(s.obj.integrated, features=top10_log2FC_markers$gene, slot='data', assay='integrated')
ggsave(filename="ureter-tumor-log2fc-heatmap.pdf", width=16,height=10,units='in')
write.csv(top_tumors_markers, file="ureter-topintVisiumMarkers.csv")

#Re-assign
sample <- "Donor4_ureter"
mb.bs <- readRDS(paste0("~/tingalab/Matt/ureter/visium/BayesSpace/dataObjects/", sample, ".Rds"))#[,read.csv(paste0("~/tingalab/Matt/human-bladder/tumors/visium/filtered-cells/", strsplit(sample, split="_")[[1]][2], ".csv"))$Barcode]
mb.int <- subset(s.obj.integrated, orig.ident==sample)
mb.bs$spatial.cluster[match(sapply(names(mb.int$seurat_clusters), FUN = function(x) strsplit(x, split="_")) %>% map_chr(c(1)), mb.bs$spot)] <- mb.int$seurat_clusters
#mb.bs <- readRDS(paste0("~/dataObjects/", sample, ".Rds"))
mb.bs %>% 
  saveRDS(file=paste0("~/ureter-dataObjects/", sample, ".Rds"))


# Enhanced Resolution Facet Highlight
# GGplot default Color Palette Emulator
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

# For 8 clusters
#pal <- gg_color_hue(8)

clusterPlotFacet <- function(sobj, cluster, palette=gg_color_hue(max(sobj$spatial.cluster))){
  pal <- rep("Grey", length(palette))
  pal[cluster] <- palette[cluster]
  clusterPlot(sobj, palette=pal)
}

scRNA <- readRDS("~/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/2020_09_21_ureter10/2020_09_21_ureter10_clustered.rds")

uro <- readRDS("~/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Uro subset/2021_04_01_ureter10_uro_PC50_res0.2/2021_04_01_ureter10_uro_PC50_res0.2_clustered.rds")
immune <- readRDS("~/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Immune subset/2021_02_25_ureter10_immune/2021_02_25_ureter10_immune_clustered.rds")
stromal <- readRDS("~/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/Stromal subset/2021_02_19_ureter10_stromal/2021_02_19_ureter10_stromal_clustered.rds")


global.anchor.labels=c("Basal Cells", "Intermediate", "Cytoxic T Cells", "T Cells", 
                       "Leukocytes", "Urothelial", "Fibroblasts", "Umbrella Cells", "Basal Cells", 
                       "Leukocytes", "NK Cells", "Endothelial Cells", "B Cells", 
                       "Smooth Muscle", "Urothelial", "Basal Cells", "Leukocytes", 
                       "Mast Cells", "Urothelial", "Fibroblasts")

global.anchor.labels.distinct=c("Basal Cells 1", "Intermediate", "Cytoxic T Cells", "T Cells", "Leukocytes 1", 
                                "Urothelial 1", "Fibroblasts 1", "Umbrella Cells", "Basal Cells 2", 
                                "Leukocytes 2", "NK Cells", "Endothelial Cells", "B Cells", 
                                "Smooth Muscle", "Urothelial 2", "Basal Cells 3", "Leukocytes 3", 
                                "Mast Cells", "Urothelial 3", "Fibroblasts 2")

urothelial.labels=c("Intermediate (CDHR5+)", "Intermediate (CRABP2-Hi)", "Basal Cells (SHH-Hi)", "Basal Cells (JUN-FOS-Hi)", 
                    "Umbrella Cells (UPK-Hi)", "Intermediate (FOXAI-Hi)", "Basal Cells (MKI67-Hi)", "Umbrella Cells (FTL-Hi)")

stromal.labels<-c("Fibroblasts (HAS1-Hi)", "Venous Endothelial Cells", "Smooth Muscle", "Fibroblasts (APOE-Hi)", "Fibroblasts (GAS1-Hi)",
                  "Fibroblasts (COL1A1-Hi)", "Arterial Endothelial Cells")
immune.labels<-c("Classical Monocytes", "Double Negative T Cells", "CD8 Trm(INFG-)", "NK Cells", "CD4 Tcm", "CD8 effector T Cells(GZMK+)", "B Cells", 
                 "Dendritic Cells", "Macrophages", "CD8 effector T Cells(GZMH+)", "Intermediate Monocytes", "Regulatory T Cells", "Non-Classical Monocytes", 
                 "Mast Cells", "CD8 Trm (INFG+)", "MAIT Cells")

anchorMapping <- function(reference, query, query.dims=15, anchor.labels,save.loc=FALSE){
  DefaultAssay(query) <- "Spatial"
  #reference@meta.data$subclass = sapply(reference@meta.data$seurat_clusters, function(i){anchor.labels[as.numeric(i)][1]})
  anchors = FindTransferAnchors(reference, query = query,  normalization.method = "LogNormalize")
  predictions.assay <- TransferData(anchorset = anchors, refdata = reference$subclass, prediction.assay=TRUE, weight.reduction=query[["pca"]], dims=1:query.dims)
  non.mapping <- c()
  for(i in 1:dim(predictions.assay)[1]){ if(sum(predictions.assay@data[i,])==0) non.mapping <- c(non.mapping, rownames(predictions.assay[i]))}
  #print("Non-integrated features: ", non.mapping)
  predictions.assay@misc$non.mapping <- non.mapping
  predictions.assay@misc$mapping <- setdiff(anchor.labels, non.mapping)
  query[["predictions"]] <- predictions.assay
  DefaultAssay(query) <- 'predictions'
  return(query)
}

s.obj.integrated <- readRDS("~/integrated-ureter-visium.Rds")
s.obj.integrated <- anchorMapping(scRNA, s.obj.integrated, query.dims=30, anchor.labels=levels(as.factor(scRNA$subclass)))
bl197129u <- subset(s.obj.integrated, orig.ident=="BL_197129U")

plot.labels <- c("Umbrella Cells","Intermediate", "Basal Cells", "Urothelial", 
                                      "Leukocytes", "Mast Cells", "Cytoxic T Cells", "T Cells", "NK Cells","B Cells",
                                      "Endothelial Cells","Smooth Muscle", "Fibroblasts"
)

scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(x) SpatialFeaturePlot(s.obj.integrated, x, images="slice1.2"))
patchwork::wrap_plots(scplots, ncol=4)


wrap_plots(p_list ,guides = 'collect', design = layout1)
# Print in 75" x 60"