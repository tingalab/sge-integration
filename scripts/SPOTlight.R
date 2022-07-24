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

ggsave(plot =  patchwork::wrap_plots(purrr::map(levels(as.factor(scRNA$subclass)), function(x) SPOTlight::spatial_scatterpie(se_obj = U1.SL,
                                                                                        cell_types_all = levels(as.factor(scRNA$subclass)),
                                                                                        cell_types_interest=x,
                                                                                        img_path = "data/U1/spatial/tissue_hires_image.png",
                                                                                        pie_scale = 0.4)), ncol = 6), filename="figures/SPOTlight/Figure_3a.pdf", device="pdf")
ggsave(plot =  patchwork::wrap_plots(purrr::map(levels(as.factor(scRNA$subclass)), function(x) SPOTlight::spatial_scatterpie(se_obj = U2.SL,
                                                                                                                             cell_types_all = cell_types_all,
                                                                                                                             cell_types_interest=x, 
                                                                                                                             img_path = "data/U2/spatial/tissue_hires_image.png",
                                                                                                                             pie_scale = 0.4), ncol = 2), file="figures/SPOTlight/Figure_3b"), device="pdf")


