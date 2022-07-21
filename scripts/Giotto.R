library(Giotto)
library(msigdbr)
library(ggplot2)
library(viridis)

results_folder = '~/Giotto/BL_190222U'

python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment(packages_to_install = c("pandas==1.1.5", "networkx==2.6.3", "python-igraph==0.9.6", "leidenalg==0.8.7",
                                                   "python-louvain==0.15", "scikit-learn==0.24.2"), force_environment = TRUE)
}

instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE)

## provide path to visium folder
data_path = '~/tingalab/Matt/ureter/visium/sr-results/BL_197129U/outs/'
data_path = '~/tingalab/Matt/ureter/visium/sr-results/BL_190222U/outs/'

x <- createGiottoVisiumObject(visium_dir = data_path, expr_data = 'raw',
                              png_name = 'tissue_hires_image.png',
                              gene_column_index = 2, instructions = instrs)


x <- updateGiottoImage(x, image_name="image", xmax_adj = 1750, ymin_adj = 1600, ymax_adj = 1750, xmin_adj = 1650) # BL_197129U
x <- updateGiottoImage(x, image_name="image", xmax_adj = 2000, ymin_adj = 1650, ymax_adj = 1650, xmin_adj = 1550) # BL_190222U

spatPlot(gobject = x, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = '2_a_spatplot_image'))

spatPlot(gobject = x, cell_color = 'in_tissue', point_size = 2,
         cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
         save_param = list(save_name = '2_c_in_tissue'))

metadata = pDataDT(x)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
x = subsetGiotto(x, cell_ids = in_tissue_barcodes)

# x <- filterGiotto(gobject = x,
#                              expression_threshold = 1,
#                              gene_det_in_min_cells = 50,
#                              min_det_genes_per_cell = 1000,
#                              expression_values = c('raw'),
#                              verbose = T)

x <- normalizeGiotto(gobject = x, scalefactor = 10000, verbose = T)
x <- addStatistics(x)
spatPlot2D(gobject = x, show_image = T, point_alpha = 0.7,
           save_param = list(save_name = '2_d_spatial_locations'))

spatPlot2D(gobject = x, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_genes', color_as_factor = F,
           save_param = list(save_name = '2_e_nr_genes'))

x <- calculateHVG(gobject = x, save_param = list(save_name = '3_a_HVGplot'))

gene_metadata = fDataDT(x)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

x <- runPCA(gobject = x, 
            genes_to_use = featgenes, 
            scale_unit = F, center = T, 
            method="factominer")

screePlot(x, ncp = 30, save_param = list(save_name = '3_b_screeplot'))

plotPCA(gobject = x,
        save_param = list(save_name = '3_c_PCA_reduction'))

## run UMAP on PCA space (default)
x <- runUMAP(x, dimensions_to_use = 1:10)
plotUMAP(gobject = x, save_param = list(save_name = '3_d_UMAP_reduction'))

## sNN network (default)
x <- createNearestNetwork(gobject = x, dimensions_to_use = 1:10, k = 15)
## Leiden clustering
x <- doLeidenCluster(gobject = x, resolution = 0.3, n_iterations = 1000)
plotUMAP(gobject = x,
         cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5,
         save_param = list(save_name = '4_a_UMAP_leiden'))

# expression and spatial
spatDimPlot(gobject = x, cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5,
            save_param = list(save_name = '5_a_covis_leiden'))

spatDimPlot(gobject = x, cell_color = 'nr_genes', color_as_factor = F,
            dim_point_size = 2, spat_point_size = 2.5,
            save_param = list(save_name = '5_b_nr_genes'))

# gini_markers_subclusters = findMarkers_one_vs_all(gobject = x,
#                                                   method = 'gini',
#                                                   expression_values = 'normalized',
#                                                   cluster_column = 'leiden_clus',
#                                                   min_genes = 20,
#                                                   min_expr_gini_score = 0.5,
#                                                   min_det_gini_score = 0.5)
# topgenes_gini = gini_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes

# violinPlot(x, genes = unique(topgenes_gini), cluster_column = 'leiden_clus',
#            strip_text = 8, strip_position = 'right',
#            save_param = list(save_name = '6_a_violinplot_gini', base_width = 5, base_height = 10))

scran_markers_subclusters = findMarkers_one_vs_all(gobject = x,
                                                   method = 'scran',
                                                   expression_values = 'normalized',
                                                   cluster_column = 'leiden_clus')
topgenes_scran = scran_markers_subclusters[, head(.SD, 2), by = 'cluster']$genes

violinPlot(x, genes = unique(topgenes_scran), cluster_column = 'leiden_clus',
           strip_text = 10, strip_position = 'right',
           save_param = list(save_name = '6_d_violinplot_scran', base_width = 5))

# cluster heatmap
plotMetaDataHeatmap(x, selected_genes = topgenes_scran,
                    metadata_cols = c('leiden_clus'),
                    save_param = list(save_name = '6_e_metaheatmap_scran'))

# umap plots
dimGenePlot(x, expression_values='normalized',
            genes = scran_markers_subclusters[, head(.SD, 1), by = 'cluster']$genes,
            cow_n_col = 3, point_size = 1,
            save_param = list(save_name = '6_f_scran_umap', base_width = 8, base_height = 5))


x <- createSpatialGrid(gobject = x,
                       sdimx_stepsize = 400,
                       sdimy_stepsize = 400,
                       minimum_padding = 0)
spatPlot(x, cell_color = 'leiden_clus', show_grid = T,
         grid_color = 'red', spatial_grid_name = 'spatial_grid', 
         save_param = list(save_name = '8_grid'))

x <- createSpatialNetwork(gobject = x, 
                          method = 'kNN', k = 5, 
                          maximum_distance_knn = 400,
                          minimum_k = 1,
                          name = 'spatial_network')

showNetworks(x)

spatPlot(gobject = x, show_network = T,
         network_color = 'blue', spatial_network_name = 'spatial_network',
         save_param = list(save_name = '9_a_knn_network'))


## rank binarization
ranktest = binSpect(x, hub_min_int = 5, bin_method="rank", spatial_network_name = "spatial_network")
spatGenePlot(x, genes = ranktest$genes[1:6], cow_n_col = 2, point_size = 1.5,
             save_param = list(save_name = '10_a_spatial_genes_km'))


## spatially correlated genes ##
ext_spatial_genes = ranktest[1:2000]$genes
ext_spatial_genes = VariableFeatures(scRNA)
cpdb1 <- read.csv("~/ureter-aux-files/cpdb_1.csv", header = FALSE) %>% unique() %>% c() %>% unlist()
cpdb2 <- read.csv("~/ureter-aux-files/cpdb_2.csv", header = FALSE) %>% unique() %>% c() %>% unlist()
cpdb <- read.csv("~/ureter-aux-files/uro_cpdb.csv", sep = "\t", header = FALSE, col.names = c("gene_a", "gene_b"))
ext_spatial_genes <- c(cpdb1, cpdb2, "SHH") %>% unique()

cpdb3 <- read.csv("~/ureter-aux-files/fibro_cpdb.csv", header = FALSE) %>% unique() %>% c() %>% unlist()
ext_spatial_genes <- c(cpdb[,1], cpdb[,2]) %>% unique()

# 1. calculate gene spatial correlation and single-cell correlation 
# create spatial correlation object
spat_cor_netw_DT = detectSpatialCorGenes(x, 
                                         method = 'network', 
                                         spatial_network_name = 'spatial_network',
                                         subset_genes = ext_spatial_genes)

# 2. identify most similar spatially correlated genes for one gene
EPCAM_top10_genes = showSpatialCorGenes(spat_cor_netw_DT, genes = 'EPCAM', show_top_genes = 10)
FGFR2_top10_genes = showSpatialCorGenes(spat_cor_netw_DT, genes = 'FGFR2', show_top_genes = 10)
spatGenePlot(x, expression_values = 'scaled',
             genes = c('EPCAM', 'HMGCS2', "KRT8", "DHRS2"), point_size = 3,
             save_param = c(save_name = '10_d_EPCAM_correlated_genes'))

spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = 3)

heatmSpatialCorGenes(x, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                     heatmap_legend_param = list(title = NULL), show_row_names = TRUE, show_cluster_annot = TRUE)

spatCor <- function(x, gene, partners, k = 2){
  geneset <- c(gene, partners)
  spat_cor_netw_DT = detectSpatialCorGenes(x, 
                                           method = 'network', 
                                           spatial_network_name = 'spatial_network',
                                           subset_genes = geneset)
  spat_cor_netw_DT = clusterSpatialCorGenes(spat_cor_netw_DT, name = 'spat_netw_clus', k = k)
  r <- spat_cor_netw_DT$cor_DT[2,]
  r$variable %<>% as.character()
  r <- r[,-c(5:8)]
  colnames(r) <- c("Gene A", "Gene B", "Spatial_corr", "Expression_corr")
  return(r)
  # plt <- heatmSpatialCorGenes(x, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', heatmap_legend_param = list(title = NULL), show_row_names = TRUE, show_cluster_annot = TRUE)
  # write.csv(spat_cor_netw_DT$cor_DT, file = paste0(gene, "_partners_correlations_uro.csv"))
  # return(plt)
}

# for(i in 1:nrow(cpdb)){
#   gene.a <- cpdb[i,1]
#   gene.b <- cpdb[i,2]
#   png(filename = paste0("./", gene.a, "x", gene.b, "_spatial_corr.png"), width = 4, height = 4, units = 'in', res = 300)
#   spatCor(x, gene.a, gene.b)
#   dev.off()
# }

z <- cbind(sapply(seq(1,nrow(cpdb)), FUN = function(i) spatCor(x, gene=cpdb[i,1], partners = cpdb[i,2])))[1:3,] %>% t 

netw_ranks = rankSpatialCorGroups(x, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                  save_param = c(save_name = '10_f_rank_correlated_groups',
                                                 base_height = 3, base_width = 5))

# Plot by cluster
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 1:3, show_top_genes = 1)
cluster_genes_DT = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_genes = 5)
cluster_genes = cluster_genes_DT$clus
names(cluster_genes) = cluster_genes_DT$gene_ID
x = createMetagenes(x, gene_clusters = cluster_genes, name = 'cluster_metagene')

spatCellPlot(x,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = netw_ranks$clusters,
             point_size = 1.5, cow_n_col = 4, 
             save_param = c(save_name = '10_g_spat_enrichment_score_plots',
                            base_width = 13, base_height = 6))

# Plot by gene 
top_netw_spat_cluster = showSpatialCorGenes(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = 1:8, show_top_genes = 1)
first_genes = top_netw_spat_cluster[, head(.SD, 1), by = clus]$gene_ID
cluster_names = top_netw_spat_cluster[, head(.SD, 1), by = clus]$clus
names(first_genes) = cluster_names
first_genes = first_genes[as.character(netw_ranks$clusters)]

spatGenePlot(x, genes = first_genes, expression_values = 'scaled', cow_n_col = 4, midpoint = 0, point_size = 2,
             save_param = c(save_name = '10_h_spat_enrichment_score_plots_genes',
                            base_width = 11, base_height = 6))

# Spatial Enrichment - Load the single cell data and create a sign matrix
scRNA <- readRDS("grouped-subsets-2.Rds")

# Rank method signature matrix - used for scRNA integration
sc_sign_matrix <- makeSignMatrixRank(
  sc_matrix = as.matrix(scRNA@assays$integrated@data),
  sc_cluster_ids = scRNA$subclass,
  ties_method = c("random"),
  gobject = NULL
)

# Using just highly variable genes
# sc_sign_matrix <- makeSignMatrixRank(
#   sc_matrix = as.matrix(scRNA@assays$RNA@data[VariableFeatqures(scRNA),]),
#   sc_cluster_ids = scRNA$subclass,
#   ties_method = c("random"),
#   gobject = NULL
# )

# Run spatial enrichment
x <- runSpatialEnrich(
  x,
  enrich_method = c("rank"),
  sign_matrix = sc_sign_matrix,
  expression_values = c("normalized"),
)


# Plot single cell enrichment patterns
spatPlot(gobject = x, cell_color=x@spatial_enrichment$rank$``, point_size = 4.5)
scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = x, cell_color=unlist(c(x@spatial_enrichment$rank[,..i])), point_size = 4.5) + ggtitle(i) + scale_fill_distiller(palette = "Spectral"))
patchwork::wrap_plots(scplots, ncol=4)
for(i in labels){
  ggsave(plot =  spatPlot(gobject = x, cell_color=unlist(c(x@spatial_enrichment$rank[,..i])), point_size = 4.5) + ggtitle(i),
         file=paste0(i, "_Giotto_scatterpie.pdf"), 
         device="pdf"
  )
}

s1 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/BL_190222U/outs/")
for(i in levels(as.factor(scRNA$subclass))){ s1@meta.data[,i] <- unlist(c(x@spatial_enrichment$rank[,..i]))}


#Plot single cell enrichment patterns using Seurat
s3 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/BL_197129U/outs/") 
for(i in levels(as.factor(scRNA$subclass))){ s3@meta.data[,i] <- unlist(c(x@spatial_enrichment$rank[,..i]))}
scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) SpatialFeaturePlot(s3, features=i, max.cutoff = 0.25) + ggtitle(i))
patchwork::wrap_plots(scplots, ncol=8)

# PAGE enrichment analysis

# Load KEGG data
kegg <- msigdbr(category="C2", subcategory = "CP:KEGG")
reactome <- msigdbr(category="C2", subcategory= "CP:REACTOME")
gs <- reactome
pathways <- unique(gs$gs_name)
gene_list <- sapply(pathways, function(x) subset(gs, gs_name == x)[,"gene_symbol"])

page_sign_matrix <- makeSignMatrixPAGE(sign_names = pathways, sign_list = gene_list)

x <- runPAGEEnrich(x , sign_matrix = page_sign_matrix, expression_values = "scaled")

spatPlot(gobject = x, cell_color=x@spatial_enrichment$PAGE$KEGG_WNT_SIGNALING_PATHWAY, point_size = 4.5) + scale_color_distiller(palette = "RdYlGn")

for(i in pathways){
  ggsave(plot =  spatPlot(gobject = x, cell_color=unlist(c(x@spatial_enrichment$PAGE[,..i])), color_as_factor = F, point_size = 4.5) + ggtitle(i),
         file=paste0(i, "_PAGE-enrich.pdf"), 
         device="pdf"
  )
}

#Plot single cell enrichment patterns using Seurat
s3 <- Load10X_Spatial("~/tingalab/Matt/ureter/visium/sr-results/BL_197129U/outs/") 
for(i in unique(scRNA$subclass)){ s3@meta.data[,i] <- unlist(c(x@spatial_enrichment$rank[,..i]))}
scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) SpatialFeaturePlot(s3, features=i) + ggtitle(i))
patchwork::wrap_plots(scplots, ncol=8)

# Cluster-cluster interaction 

CPScore <-cellProximityEnrichment(
  gobject=x,
  spatial_network_name = "spatial_network",
  cluster_column="leiden_clus",
  number_of_simulations = 1000,
  adjust_method = "bonferroni"
)

cellProximityBarplot(
  gobject=x,
  CPscore=CPScore,
  min_orig_ints = 5,
  min_sim_ints = 5,
  p_val = 0.05,
  save_plot = NA,
  save_param = list(),
  default_save_name = "cellProximityBarplot"
)

x <- cellProximityNetwork(x, CPScore = CPScore,  remove_self_edges = T, only_show_enrichment_edges = T)

# Carry over Seurat object data into Giotto format

y <- createGiottoObject(raw_exprs = s3[["Spatial"]]@counts, norm_expr = s3[["Spatial"]]@data, norm_scaled_expr = s3[["Spatial"]]@scale.data, 
)


x <- readRDS("~/Giotto/Giotto-w-seurat-counts.Rds")
x <- readRDS(file = "~/Giotto/cell-type-annotated-Giotto.Rds")