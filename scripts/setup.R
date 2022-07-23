library(Seurat)
library(dplyr)
library(tidyverse)
library(BayesSpace)
library(purrr)
library(janitor)
library(magrittr)
library(purrr)
library(patchwork)
library(stringr)

# Set this to your sge-integration folder
setwd("/home/bradlem4/sge-integration/")

# Download Visium data from GitHub release
download.file(url = "https://github.com/tingalab/sge-integration/releases/download/V1/protocol-data.tar.gz", destfile = "protocol-data.tar.gz")
untar("data.tar.gz")
file.remove("data.tar.gz")

# Download single-cell data from GEO
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184111/suppl/GSE184111_scRNA_ureter10_normalized_counts.txt.gz", destfile = "data/scRNA/count-matrix.txt.gz")
system("gunzip data/scRNA/count-matrix.txt.gz")
scRNA <- read.csv(file = "data/scRNA/count-matrix.txt", sep = "\t")
colnames(scRNA) %<>% str_replace_all(., pattern = "\\.", replacement="-")

# Load metadata, remove some cells that don't have associated metadata. Then, create a Seurat object.
md <- read.csv(file = "data/scRNA/scRNA-metadata.csv", row.names = 1)
scRNA <- scRNA[,rownames(md)] 
scRNA %<>% CreateSeuratObject(meta.data = md)
#scRNA@meta.data %<>% cbind(.,md)

# Save the Seurat object for use with other integration methods.
saveRDS(scRNA, file = "data/scRNA/ureter-scRNA.Rds")