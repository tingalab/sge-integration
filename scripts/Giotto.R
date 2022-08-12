library(Giotto)
library(msigdbr)
library(ggplot2)
library(viridis)

# Set this to your sge-integration folder
setwd("/home/bradlem4/sge-integration/")

source("scripts/functions.R")
scRNA <- readRDS("data/scRNA/ureter-scRNA.Rds")

#----- Configure workspace with Giotto

results_folder = 'figures/Giotto'

python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment(packages_to_install = c("pandas==1.1.5", "networkx==2.6.3", "python-igraph==0.9.6", "leidenalg==0.8.7",
                                                   "python-louvain==0.15", "scikit-learn==0.24.2"), force_environment = TRUE)
}
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE)

## Create new Giotto objects for each set of Visium data. We have performed the adjustment necessary, however new sets of Visium data will need to adjust the min_adj and max_adj values as necessary.
# For more details, please see the Giotto vignette. 

U1 <- createGiottoVisiumObject(visium_dir = 'data/U1', expr_data = 'filter', 
                              h5_visium_path = 'data/U1/filtered_feature_bc_matrix.h5',
                              h5_tissue_positions_path = 'data/U1/spatial/tissue_positions_list.csv', 
                              h5_image_png_path = 'data/U1/spatial/tissue_lowres_image.png',
                              h5_json_scalefactors_path= 'data/U1/spatial/scalefactors_json.json',
                              gene_column_index = 2, instructions = instrs,
                              xmax_adj = 3200, ymin_adj = 2600, ymax_adj = 2900, xmin_adj = 1600)

U2 <-  createGiottoVisiumObject(visium_dir = 'data/U2', expr_data = 'filter', 
                                     h5_visium_path = 'data/U2/filtered_feature_bc_matrix.h5',
                                     h5_tissue_positions_path = 'data/U2/spatial/tissue_positions_list.csv', 
                                     h5_image_png_path = 'data/U2/spatial/tissue_lowres_image.png',
                                     h5_json_scalefactors_path= 'data/U2/spatial/scalefactors_json.json',
                                     gene_column_index = 2, instructions = instrs, xmax_adj = 2000, ymin_adj = 1500, ymax_adj = 1600, xmin_adj = 1400
                                      )

# These functions can help you determine whether your alignment is correct. See troubleshooting section of the STAR protocol for more details. 
# spatPlot(gobject = U1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
#          save_param = list(save_name = 'U1_spatplot_image'))
# 
# spatPlot(gobject = U2, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
#          save_param = list(save_name = 'U2_spatplot_image'))

U1.Giotto<- preProcessGiotto(U1, "U1")
U2.Giotto<-preProcessGiotto(U2, "U2")

# First you have to create the rank method signature matrix.
sc_sign_matrix <- makeSignMatrixRank(
  sc_matrix = as.matrix(scRNA@assays$RNA@data),
  sc_cluster_ids = scRNA$subclass,
  ties_method = c("random"),
  gobject = NULL
)


# Run spatial enrichment for both samples
U1.Giotto <- runSpatialEnrich(
  U1.Giotto,
  enrich_method = c("rank"),
  sign_matrix = sc_sign_matrix,
  expression_values = c("normalized"),
)

U2.Giotto <- runSpatialEnrich(
  U2.Giotto,
  enrich_method = c("rank"),
  sign_matrix = sc_sign_matrix,
  expression_values = c("normalized"),
)

# Plot single cell enrichment patterns and simultaneously save figures in folder.
scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = U1.Giotto, cell_color=unlist(c(U1@spatial_enrichment$rank[,..i])), point_size = 2) + ggtitle(i) + scale_fill_distiller(palette = "Spectral"))
patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Giotto/Figure_2a.pdf", width = 25, height = 20, units = "in", dpi = 300)

scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = U2.Giotto, cell_color=unlist(c(U2@spatial_enrichment$rank[,..i])), point_size = 2) + ggtitle(i) + scale_fill_distiller(palette = "Spectral"))
patchwork::wrap_plots(scplots, ncol=4)%T>% ggsave(filename = "figures/Giotto/Figure_2b.pdf", width = 25, height = 20, units = "in", dpi = 300)
                      
# OPTIONAL: Use the following code for printing larger fonts for figures

scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = U1.Giotto, cell_color=unlist(c(U1.Giotto@spatial_enrichment $cell$rna $rank[,..i])), point_size = 2) + 
theme(title = element_text(size=18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          axis.title = element_text(size=15),
          axis.text = element_text(size=15) + 
ggtitle(i) + scale_fill_distiller(palette = “Spectral”))

patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Giotto/Figure_2a.pdf", width = 25, height = 20, units = "in", dpi = 300)

                      
scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = U2.Giotto, cell_color=unlist(c(U2.Giotto@spatial_enrichment $cell$rna $rank[,..i])), point_size = 2) + 
theme(title = element_text(size=18),
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          axis.title = element_text(size=15),
          axis.text = element_text(size=15) + 
ggtitle(i) + scale_fill_distiller(palette = “Spectral”))

patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Giotto/Figure_2b.pdf", width = 25, height = 20, units = "in", dpi = 300)

