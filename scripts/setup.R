
library(Seurat)
library(dplyr)
library(tidyverse)
library(purrr)
library(janitor)
library(magrittr)
library(patchwork)
library(stringr)
library(R.utils)

# Set sge-integration folder as the currect working directory

# replace the following path within quotes to your current/home directory where you downloaded the sge-integration repository
homedir<-"/home/sonas/star_protocol" 
setwd(paste0(homedir,"/sge-integration/"))

# Download Visium data from GitHub release
download.file(url = "https://github.com/tingalab/sge-integration/releases/download/V1/protocol-data.tar.gz", destfile = "protocol-data.tar.gz")
gunzip("protocol-data.tar.gz")
untar("protocol-data.tar")
file.remove("protocol-data.tar")

# The above lines will download protocol data packaged within "data" named folder to your sge-integration folder #
# This "data" directory harbors subdirectories "U1" and "U2" that contain Visium spatial transcriptomics data for 2 ureter samples
# Another sudirectory "scRNA" contains single cell RNA-seq metadata from our Ureter manuscript (Fink and Sona et al., 2022)
# The matching single cell RNA-seq data can be downloaded to this "scRNA" subdirectory using the following code:

# Download single-cell data from GEO
download.file(url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184111/suppl/GSE184111_scRNA_ureter10_normalized_counts.txt.gz", destfile = "data/scRNA/count-matrix.txt.gz")
gunzip("data/scRNA/count-matrix.txt.gz")

scRNA <- read.csv(file = "data/scRNA/count-matrix.txt", sep = "\t")
colnames(scRNA) %<>% str_replace_all(., pattern = "\\.", replacement="-")

# Load metadata, remove some cells that don't have associated metadata. Then, create a Seurat object.
md <- read.csv(file = "data/scRNA/scRNA-metadata.csv", row.names = 1)
scRNA <- scRNA[,rownames(md)] 
scRNA %<>% CreateSeuratObject(meta.data = md)
#scRNA@meta.data %<>% cbind(.,md)

# Save the Seurat object for use with other integration methods.
saveRDS(scRNA, file = "data/scRNA/ureter-scRNA.Rds")
