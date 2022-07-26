  library(Seurat)
  library(gt)
  library(dplyr)
  library(Matrix)
  library(data.table)
  library(SPOTlight)
  library(igraph)
  library(RColorBrewer)
  
  # Set this to your sge-integration folder
  setwd("/home/bradlem4/sge-integration/")
  
  source("scripts/functions.R")
  scRNA <- readRDS("data/scRNA/ureter-scRNA.Rds")

#----- Load data

# Since SPOTlight uses the same format as Seurat, we can repurpose the methods to load
U1.SL <- preProcessSeuratVisium("data/U1", normalization = "LogNormalize")
U2.SL <- preProcessSeuratVisium("data/U2", normalization = "LogNormalize")

# You may need to generate your own markers. This will take much longer if you do we give the option of loading some from a CSV file for the sample data.
#markers <- read.csv("data/scRNA/markers.csv")
Idents(scRNA) <- scRNA$subclass
markers <- FindAllMarkers(scRNA, group.by = "subclass")

U1.SL <- spotlightDeconvolve(U1.SL, scRNA, markers)
U2.SL <- spotlightDeconvolve(U2.SL, scRNA, markers)

cell.names <- names(which(colSums(U1.SL@meta.data[,-c(1:4, ncol(U1.SL@meta.data))]) > 0))
ggsave(plot =  patchwork::wrap_plots(purrr::map(cell.names, function(x) SPOTlight::spatial_scatterpie(se_obj = U1.SL,
    cell_types_all = cell.names,
    cell_types_interest=x,
    img_path = "data/U1/spatial/tissue_lowres_image.png",
    pie_scale = 0.4)), ncol = 4), filename="figures/SPOTlight/Figure_3a.pdf", device="pdf", width = 25, height = 25, units = "in", dpi = 300)

cell.names <- names(which(colSums(U2.SL@meta.data[,-c(1:4, ncol(U2.SL@meta.data))]) > 0))
ggsave(plot =  patchwork::wrap_plots(purrr::map(cell.names, function(x) SPOTlight::spatial_scatterpie(se_obj = U2.SL,
    cell_types_all = cell.names,
     cell_types_interest=x, 
     img_path = "data/U2/spatial/tissue_lowres_image.png",
     pie_scale = 0.4)), ncol = 4), filename="figures/SPOTlight/Figure_3b.pdf", 
    device="pdf", width = 25, height = 25, units = "in", dpi = 300)


