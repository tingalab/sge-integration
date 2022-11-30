library(Matrix)
library(data.table)
library(SPOTlight)
library(igraph)
library(RColorBrewer)
  
# Set this to your sge-integration folder
mydir='/home/sonas/star_protocol/'

setwd(paste0(mydir,"sge-integration/"))
  
source("scripts/functions.R")
scRNA <- readRDS("data/scRNA/ureter-scRNA.Rds")

#----- Load data

# Since SPOTlight uses the same format as Seurat, we can repurpose the methods to load
U2.SL <- preProcessSeuratVisium("data/U2", normalization = "LogNormalize")

# You may need to generate your own markers. This will take much longer if you do we give the option of loading some from a CSV file for the sample data.
#markers <- read.csv("data/scRNA/markers.csv")
Idents(scRNA) <- scRNA$subclass
markers <- FindAllMarkers(scRNA, group.by = "subclass")

#U2.SL <- spotlightDeconvolve(U2.SL, scRNA, markers)

genes<-read.csv(file ="data/scRNA/genes.csv")[,1]
U2.SL <- spotlightDeconvolve(vis = U2.SL, scrna = scRNA,markers =  markers,hvg=genes)


cell.names <- names(which(colSums(U2.SL@meta.data[,-c(1:4, ncol(U2.SL@meta.data))]) > 0))
ggsave(plot =  patchwork::wrap_plots(purrr::map(cell.names, function(x) SPOTlight::spatial_scatterpie(se_obj = U2.SL,
    cell_types_all = cell.names,
     cell_types_interest=x, 
     img_path = "data/U2/spatial/tissue_lowres_image.png",
     pie_scale = 0.4)), ncol = 4), filename="figures/SPOTlight/Figure_3b.pdf", 
    device="pdf", width = 25, height = 25, units = "in", dpi = 300)

# OPTIONAL: Added another version of plotting function to print larger fonts in the figures
                                     

cell.names <- names(which(colSums(U2.SL@meta.data[,-c(1:4, ncol(U2.SL@meta.data))]) > 0))
scplots<-purrr::map(cell.names, function(x) SPOTlight::spatial_scatterpie(se_obj = U2.SL,  cell_types_all = cell.names,  cell_types_interest=x,  
img_path = "data/U2/spatial/tissue_lowres_image.png",
pie_scale = 0.4)+
 theme(title = element_text(size=20),
       legend.text = element_text(size = 18),
       legend.title = element_text(size = 18),
       axis.title = element_text(size=18),
       axis.title.y = element_text(angle = 90),
       axis.text = element_text(size=18)))

patchwork::wrap_plots(scplots, ncol=3) %T>% ggsave(filename="figures/SPOTlight/Figure3b.pdf", device="pdf", width = 25, height = 30, units = "in", dpi = 300)

