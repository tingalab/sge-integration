library(Seurat)
library(gt)
library(dplyr)
library(Matrix)
library(data.table)
library(SPOTlight)
library(igraph)
library(RColorBrewer)

loadVisium <- function(path){
  sample <- Load10X_Spatial(path) #Load10X_Spatial("~/tingalab/Matt/visium_results/BL_197129U/outs/")
  sample <- PercentageFeatureSet(sample, pattern="^MT-", col.name="percent.mt")
  sample <- SCTransform(sample, assay="Spatial", vars.to.regress = "percent.mt", return.only.var.genes = FALSE) %>%
    RunPCA()
  return(sample)
}

BL1971 <- loadVisium("~/tingalab/Matt/ureter/visium/sr-results/BL_197129U/outs/")
scRNA <- readRDS("~/tingalab/Manuscripts/2021_UreterManuscript/ANALYSIS/Seurat/2020_09_21_ureter10/2020_09_21_ureter10_clustered.rds") %>%
  SCTransform(vars.to.regress = "perc.mt")


ElbowPlot(BL1971)
UMAP <- function(sample, n_dims=25, res=0.3){
  sample <- FindNeighbors(sample,reduction="pca", dims=1:n_dims)
  sample <- FindClusters(sample, resolution=res)
  sample <- RunUMAP(sample,reduction="pca", dims=1:n_dims)
  return(sample)
}

global.anchor.labels.distinct=c("Basal Cells 1", "Intermediate", "Cytoxic T Cells", "T Cells", "Leukocytes 1", 
                                "Urothelial 1", "Fibroblasts 1", "Umbrella Cells", "Basal Cells 2", 
                                "Leukocytes 2", "NK Cells", "Endothelial Cells", "B Cells", 
                                "Smooth Muscle", "Urothelial 2", "Basal Cells 3", "Leukocytes 3", 
                                "Mast Cells", "Urothelial 3", "Fibroblasts 2")

# global.anchor.labels=c("Basal Cells", "Intermediate", "Cytoxic T Cells", "T Cells", 
#                        "Leukocytes", "Urothelial", "Fibroblasts", "Umbrella Cells", "Basal Cells", 
#                        "Leukocytes", "NK Cells", "Endothelial Cells", "B Cells", 
#                        "Smooth Muscle", "Urothelial", "Basal Cells", "Leukocytes", 
#                        "Mast Cells", "Urothelial", "Fibroblasts")

global.anchor.labels.reclustered=c("T Cells", "Intermediate (Krt13+)", "Leukocytes 1", "Cytoxic T Cells", "Basal Cells 1", "Fibroblasts 1",
                                   "Umbrella Cells", "Urothelial 1", "Intermediate (Krt8+)", "Intermediate (Jun/Fos-Hi, Krt13+)", "Leukocytes 2",
                                   "NK Cells", "Endothelial Cells", "B Cells", "Smooth Muscle Cells", "Urothelial 2", "Leukocytes 3", "Basal Cells 3", "Mast Cells", "Fibroblasts 2")

anchor.labels <- global.anchor.labels.reclustered
scRNA@meta.data$new_subclass = sapply(scRNA@meta.data$new_clusters, function(i){anchor.labels[as.numeric(i)][1]})


scRNA@meta.data %>%
  dplyr::count(subclass) %>%
  gt::gt(.[-1, ]) %>%
  gt::tab_header(
    title = "Cell types present in the reference dataset",
  ) %>%
  gt::cols_label(
    subclass = gt::html("Cell Type")
  )

Seurat::Idents(object = scRNA) <- scRNA@meta.data$subclass
cluster_markers_all <- Seurat::FindAllMarkers(object = scRNA, 
                                              slot = "data",
                                              verbose = TRUE, 
                                              only.pos = TRUE)
saveRDS(cluster_markers_all, file = "scRNA_cluster_markers-new-clusters.Rds")
write.csv(cluster_markers_all, file="scRNA_cluster_markers-reclustered.csv")

spotlight_ls <- spotlight_deconvolution(
  se_sc = scRNA,
  counts_spatial = s.obj.integrated@assays$Spatial@counts,
  clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
  cluster_markers = cluster_markers_all, # Dataframe with the marker genes
  cl_n = 102, # number of cells per cell type to use
  hvg = 4273, # Number of HVG to use
  ntop = NULL, # How many of the marker genes to use (by default all)
  transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
  method = "nsNMF", # Factorization method
  min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
)

nmf_mod <- spotlight_ls[[1]]
decon_mtrx <- spotlight_ls[[2]]

decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
rownames(decon_mtrx) <- colnames(s.obj.integrated)

decon_df <- decon_mtrx %>%
  data.frame() %>%
  tibble::rownames_to_column("barcodes")

s.obj.integrated@meta.data <- s.obj.integrated@meta.data %>%
  tibble::rownames_to_column("barcodes") %>%
  dplyr::left_join(decon_df, by = "barcodes") %>%
  tibble::column_to_rownames("barcodes")

cell_types_all <- head(names(which(colSums(decon_mtrx) > 0)), -1)

SPOTlight::spatial_scatterpie(se_obj = s.obj.integrated,
                              cell_types_all = cell_types_all,
                               cell_types_interest="Endothelial.Cells",
                              slice="slice1.2",
                              img_path = "~/tingalab/Matt/ureter/visium/sr-results/BL_197129U/outs/spatial/tissue_lowres_image.png",
                              pie_scale = 0.4) + ggtitle(label="Endothelial Cells")

sample.names <- Images(s.obj.integrated)

for(i in cell_types_all){
  ggsave(plot =  patchwork::wrap_plots(purrr::map(sample.names, function(x) SPOTlight::spatial_scatterpie(se_obj = s.obj.integrated,
                                                                                        cell_types_all = cell_types_all,
                                                                                        cell_types_interest=i, 
                                                                                        slice=x,
                                                                                        img_path = paste0("~/tingalab/Matt/ureter/visium/sr-results/", x, "/outs/spatial/tissue_lowres_image.png"),
                                                                                        pie_scale = 0.4) + ggtitle(label=x)), ncol = 2),
         file=paste0(i, "_scatterpie.pdf"), 
         device="pdf"
  )
}

                              
patchwork::wrap_plots(spotlight.plots, ncol=5)


# Spatial interaction graph

graph_ntw <- SPOTlight::get_spatial_interaction_graph(decon_mtrx = decon_mtrx[, cell_types_all])


deg <- degree(graph_ntw, mode="all")

# Get color palette for difusion
edge_importance <- E(graph_ntw)$importance

# Select a continuous palette
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'seq',]

# Create a color palette
getPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))

# Get how many values we need
grad_edge <- seq(0, max(edge_importance), 0.1)

# Generate extended gradient palette dataframe
graph_col_df <- data.frame(value = as.character(grad_edge),
                           color = getPalette(length(grad_edge)),
                           stringsAsFactors = FALSE)
# Assign color to each edge
color_edge <- data.frame(value = as.character(round(edge_importance, 1)), stringsAsFactors = FALSE) %>%
  dplyr::left_join(graph_col_df, by = "value") %>%
  dplyr::pull(color)

# Open a pdf file
plot(graph_ntw,
     # Size of the edge
     edge.width = edge_importance,
     edge.color = color_edge,
     # Size of the buble
     vertex.size = deg,
     vertex.color = "#cde394",
     vertex.frame.color = "white",
     vertex.label.color = "black",
     vertex.label.family = "Ubuntu", # Font family of the label (e.g.“Times”, “Helvetica”)
     layout = layout.circle)