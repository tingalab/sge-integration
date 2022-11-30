library(Seurat)
library(dplyr)
library(tidyverse)
library(purrr)
library(janitor)
library(magrittr)
library(patchwork)
library(stringr)

#--------Setup-------

# This should be your directory/path where your sge-integration repository resides
mydir='/home/sonas/star_protocol/'

# Set this to your sge-integration folder
setwd(paste0(mydir,"sge-integration/"))

# Source the functions file to make them available to our active session of R
source("scripts/functions.R")


# Save the Seurat object for use with other integration methods.
scRNA <- readRDS(file = "data/scRNA/ureter-scRNA.Rds")
genes<-read.csv(file ="data/scRNA/genes.csv")[,1]


## ----- Load Visium samples ------
# This block of code will do a few things, but in short it will prepare the individual Visium samples for integration.
# If your single-cell data is normalized using SCTransform then you need to normalize your Visium data using SCTransform. In short, they just have to match. In order to SCT normalize your data, put "SCT" instead of "LogNormalize" below. 

U2.Seurat <- preProcessSeuratVisium("data/U2", normalization = "LogNormalize")

#----- Perform anchor mapping----

# This function will find "anchor" genes that are represented prominently in both the scRNA and SGE data. It will then use these to help approximate the distribution of single-cell identities in the spatial data.
U2.Seurat <- anchorMapping(scRNA, U2.Seurat, feats = genes, query.dims=30, anchor.labels=levels(as.factor(scRNA$subclass)))

# Plot the predictions made by the single-cell data

scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(x) SpatialFeaturePlot(U2.Seurat, x))
patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Seurat/Figure_1b.pdf", width =25, height = 25, units = "in", dpi = 300)
     
# OPTIONAL: To generate plots with larger fonts use the following set of code instead

scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(x) SpatialFeaturePlot(U2.Seurat, x) +
                      theme(legend.key.size = unit(10, "mm"),
                            legend.text = element_text(size = 15),
                            legend.title = element_text(size = 20)))
patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Seurat/Figure_1b.pdf", width =25, height = 25, units = "in", dpi = 300)

                      
                    
  

