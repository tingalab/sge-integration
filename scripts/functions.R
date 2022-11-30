library(Seurat)
library(dplyr)
library(tidyverse)
library(janitor)
library(magrittr)
library(purrr)
library(patchwork)

# This file contains all the pre-written functions for SGE integration with single-cell data. 
# This file should be sourced with a command at the beginning of the protocol.

preProcessSeuratVisium<- function(file, normalization = "SCT"){
    if(normalization == "SCT"){
      Load10X_Spatial(file) %>% 
        PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
        SCTransform(assay="Spatial", vars.to.regress = "percent.mt", return.only.var.genes = FALSE) %>%
        RunPCA() %>%
        return()
    }else if(normalization == "LogNormalize"){
      Load10X_Spatial(file) %>% 
        PercentageFeatureSet(pattern="^MT-", col.name="percent.mt") %>%
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
        FindVariableFeatures(nfeatures = 3000) %>%
        ScaleData(vars.to.regress = "percent.mt") %>%
        RunPCA() %>%
        return()
    } else {
      warning("Please pick either SCT or LogNormalize for normalization.")
    }
}

UMAP <- function(sample, n_dims=15, res=0.5){
  sample <- FindNeighbors(sample,reduction="pca", dims=1:n_dims)
  sample <- FindClusters(sample, resolution=res)
  sample <- RunUMAP(sample,reduction="pca", dims=1:n_dims)
  return(sample)
}

anchorMapping <- function(reference, query, feats, query.dims=15, normalization.method = "LogNormalize", anchor.labels,save.loc=FALSE){
  DefaultAssay(query) <- "Spatial"
  #reference@meta.data$subclass = sapply(reference@meta.data$seurat_clusters, function(i){anchor.labels[as.numeric(i)][1]})
  anchors = FindTransferAnchors(reference, query = query,  normalization.method = normalization.method, features = feats)
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

#------ Giotto -------

preProcessGiotto<-function(gobject, name){
  metadata = pDataDT(gobject)
  in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
  gobject <- subsetGiotto(gobject, cell_ids = in_tissue_barcodes)
  gobject <- normalizeGiotto(gobject, scalefactor = 10000, verbose = T)
  gobject <- addStatistics(gobject)
  gobject <- calculateHVG(gobject, save_param = list(save_name = paste0(name, '_HVGplot')))
  gene_metadata = fDataDT(gobject)
  featgenes = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID
  gobject <- runPCA(gobject, 
              genes_to_use = featgenes, 
              scale_unit = F, center = T, 
              method="factominer")
  
  gobject <- createSpatialGrid(gobject,
                         sdimx_stepsize = 400,
                         sdimy_stepsize = 400,
                         minimum_padding = 0)
  
  gobject <- createSpatialNetwork(gobject = gobject, 
                            method = 'kNN', k = 5, 
                            maximum_distance_knn = 400,
                            minimum_k = 1,
                            name = paste0(name, 'spatial_network'))
  return(gobject)
}

# ---------- SPOTlight -------

# spotlightDeconvolve <- function(vis, scrna, markers, cellsper=107,hvg=4273){
#   spotlight_ls <- spotlight_deconvolution(
#     se_sc = scrna,
#     counts_spatial = vis@assays$Spatial@counts,
#     clust_vr = "subclass", # Variable in sc_seu containing the cell-type annotation
#     cluster_markers = markers, # Dataframe with the marker genes
#     cl_n = cellsper, # number of cells per cell type to use
#     hvg = hvg, # Number of HVG to use
#     transf = "uv", # Perform unit-variance scaling per cell and spot prior to factorzation and NLS
#     method = "nsNMF", # Factorization method
#     min_cont = 0 # Remove those cells contributing to a spot below a certain threshold 
#   )
#   
#   nmf_mod <- spotlight_ls[[1]]
#   decon_mtrx <- spotlight_ls[[2]]
#   
#   decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
#   decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
#   decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
#   rownames(decon_mtrx) <- colnames(vis)
#   
#   decon_df <- decon_mtrx %>%
#     data.frame() %>%
#     tibble::rownames_to_column("barcodes")
#   
#   vis@meta.data <- vis@meta.data %>%
#     tibble::rownames_to_column("barcodes") %>%
#     dplyr::left_join(decon_df, by = "barcodes") %>%
#     tibble::column_to_rownames("barcodes")
#   return(vis) 
# }


spotlightDeconvolve <- function(vis, scrna, markers,hvg){
  spotlight_ls <- SPOTlight(
    x = scrna,
    y=vis,
    mgs = markers, # Dataframe with the marker genes
    hvg = hvg, # Highly variable genes
    groups = as.character(scrna$subclass),
    model = "ns",
    gene_id='gene',
    group_id = 'cluster',
    weight_id='avg_log2FC',
    assay_sc='RNA',
    assay_sp='Spatial',
  )
  
  nmf_mod <- spotlight_ls[[1]]
  decon_mtrx <- spotlight_ls[[2]]
  
  decon_mtrx_sub <- decon_mtrx[, colnames(decon_mtrx) != "res_ss"]
  decon_mtrx_sub[decon_mtrx_sub < 0.08] <- 0
  decon_mtrx <- cbind(decon_mtrx_sub, "res_ss" = decon_mtrx[, "res_ss"])
  rownames(decon_mtrx) <- colnames(vis)
  
  decon_df <- decon_mtrx %>%
    data.frame() %>%
    tibble::rownames_to_column("barcodes")
  
  vis@meta.data <- vis@meta.data %>%
    tibble::rownames_to_column("barcodes") %>%
    dplyr::left_join(decon_df, by = "barcodes") %>%
    tibble::column_to_rownames("barcodes")
    return(vis) 
}
