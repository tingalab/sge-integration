This notebook provides all the instructions to run scRNA integration
into visium data as described in our STAR protocol **“Methods for
single-cell RNA-seq and spatial gene expression integration and
interactive visualization”**. The users can also find all the necessary
instructions in our protocol and all the necessary scripts in our github
repository <https://github.com/tingalab/sge-integration>

# Before you begin

This section provides the setup steps (installation and data download)

## 1. Download Github repository

For Mac users: Open terminal and browse to your preferred directory or
location to download the github repository.


    # Replace mydir path with your own path
    mydir='/home/sonas/star_protocol/'
    cd ${mydir}

    git clone https://github.com/tingalab/sge-integration.git

This will download the entire sge-integration github repository to your
current directory, with *scripts* folder (containing all scripts
described here), *shiny-app* folder (containing shiny app related
scripts) and *figures* folder (with *Giotto*, *SPOTlight* and *Seurat*
as blank subfolders) to save your output from the respective workflows.
You can view the repository structure using the following code:


    mydir='/home/sonas/star_protocol/'

    cd ${mydir}/sge-integration
    ls -d -- */* 

    ## data/scRNA
    ## data/U1
    ## data/U2
    ## figures/Giotto
    ## figures/Seurat
    ## figures/SPOTlight
    ## scripts/functions.R
    ## scripts/Giotto.R
    ## scripts/setup.R
    ## scripts/Seurat.R
    ## scripts/SPOTlight.R
    ## shiny-app/app.R

*Note: b. If you are on a Windows platform, download the code in the
form of a .zip file from the repository’s web page. After downloading
the zip file, decompress the file and ensure all the contents are
identical to the GitHub repository.*

## 2. Install tools and dependencies

Please use our script `install_dependencies.R` that contains the code
and instructions for installation of dependencies. We are providing the
code here as well.

    install.packages("Seurat", version="4.2.0")
    install.packages("dplyr", version="v1.0.8")
    install.packages("tidyverse", version="v1.3.1")
    install.packages("purrr", version="v0.3.4")
    install.packages("janitor", version="v2.1.0")
    install.packages("magrittr", version="v2.0.2")
    install.packages("patchwork", version="v1.1.1")
    install.packages("ggplot2", version="v3.3.6")
    install.packages("viridis", version="v0.6.2")
    install.packages("Matrix", version="v1.4-0")
    install.packages("igraph", version="v1.3.1")
    install.packages("RColorBrewer", version="v1.1-3")
    install.packages("shiny", version="v1.7.1")
    install.packages("miniUI", version="v1.1.1.1")
    install.packages("layer", version="v0.0.1")
    install.packages("grid", version="v4.0.3")
    install.packages("Cairo", version="v1.5-15")
    install.packages("R.utils", version="v2.12")
    install.packages("devtools")
    install.packages("remotes",version="2.4.2")

    # Now install other packages (Giotto and SPOTlight)

    # For SPOTlight installation and other documentation, please refer to https://marcelosua.github.io/SPOTlight/
    devtools::install_github("https://github.com/MarcElosua/SPOTlight")


    # For detailed installation instructions and troubleshooting please see https://giottosuite.readthedocs.io/en/latest/gettingstarted.html

    #remotes::install_github("drieslab/Giotto@suite")
    remotes::install_github("RubD/Giotto@v1.1.0") 

    # installGiottoEnvironment(packages_to_install = c("pandas==1.1.5", "networkx==2.6.3", "python-igraph==0.9.6", "leidenalg==0.8.7","python-louvain==0.15", "scikit-learn==0.24.2"), force_environment = TRUE)
    #            
    installGiottoEnvironment(packages_to_install = c("pandas", 
                                                     "networkx", 
                                                     "python-igraph", 
                                                     "leidenalg",
                                                     "python-louvain",
                                                     "scikit-learn"), 
                             force_environment = TRUE)

## 3. Run Setup.R (for data download)

Please run `setup.R` script (provided in the script folder of the
`sge-integration` repository) to download test data. A stepwise
breakdown is provided here:

    # If running a new session, please run these two lines to load packages from your R library
    mylib="/home/sonas/R/mylib/" # replace this with your own R library path
    .libPaths(c(mylib, .libPaths()))

    # Load depedencies
    library(Seurat)
    library(dplyr)
    library(tidyverse)
    library(BayesSpace)
    library(purrr)
    library(janitor)
    library(magrittr)
    library(patchwork)
    library(stringr)
    library(R.utils)


    homedir<-"home/sonas/star_protocol" # replace the path within quotes to your current/home directory where you downloaded the sge-integration repository

    # optional: setwd or add full path in subsequent sections
    setwd(paste0(homedir,"/sge-integration/"))

Start downloading data

    # Download Visium data from GitHub release
    download.file(url = "https://github.com/tingalab/sge-integration/releases/download/V1/protocol-data.tar.gz", destfile = "protocol-data.tar.gz")
    gunzip("protocol-data.tar.gz")
    untar("protocol-data.tar")
    file.remove("protocol-data.tar")

The above lines will download protocol data packaged within `data` named
folder to your sge-integration folder. The `data` directory has 3
subdirectories: `scRNA`, `U1` and `U2`.

You should have the following meta data files in your `scRNA` data
folder


    # This should be your directory/path
    mydir='/home/sonas/star_protocol/'

    cd ${mydir}/sge-integration
    ls -d -- data/scRNA/*

    ## data/scRNA/count-matrix.txt
    ## data/scRNA/genes.csv
    ## data/scRNA/markers.csv
    ## data/scRNA/scRNA-metadata.csv
    ## data/scRNA/ureter-scRNA.Rds

The U1 and U2 contain visium data (anayzed by Seurat), organized as
follows:


    # This should be your directory/path
    mydir='/home/sonas/star_protocol/'

    cd ${mydir}/sge-integration
    ls -d -- data/U1/*/*

    ## data/U1/spatial/aligned_fiducials.jpg
    ## data/U1/spatial/detected_tissue_image.jpg
    ## data/U1/spatial/scalefactors_json.json
    ## data/U1/spatial/Thumbs.db
    ## data/U1/spatial/tissue_hires_image.png
    ## data/U1/spatial/tissue_lowres_image.png
    ## data/U1/spatial/tissue_positions_list.csv

Now download, process and save single cell data from GEO

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

The `scRNA` subdirectory should now contain all of the following data:



    # This should be your directory/path
    mydir='/home/sonas/star_protocol/'

    cd ${mydir}/sge-integration
    ls -d -- data/scRNA/*

    ## data/scRNA/count-matrix.txt
    ## data/scRNA/genes.csv
    ## data/scRNA/markers.csv
    ## data/scRNA/scRNA-metadata.csv
    ## data/scRNA/ureter-scRNA.Rds

# Step-by-step method details

## Perform integration in Seurat

To perform integration using Seurat, please run the script `Seurat.R`.
The code and brief description for it are provided below:

**Load scRNA and visium data**

-   Please note that for the protocol demo, we will just load and
    process just 1 Visium sample `U2`

-   If your single-cell data is normalized using SCTransform then you
    need to normalize your Visium data using SCTransform. In short, they
    just have to match. In order to SCT normalize your data, put “SCT”
    instead of “LogNormalize” below.

<!-- -->

    # If running a new session (recommended), please run these two lines to load packages from your R library
    mylib="/home/sonas/R/mylib/" # replace this with your own R library path
    .libPaths(c(mylib, .libPaths()))

    #--------Setup-------

    # This should be your directory/path where your sge-integration repository resides
    mydir='/home/sonas/star_protocol/'

    # Set this to your sge-integration folder
    setwd(paste0(mydir,"sge-integration/"))

    # Save the Seurat object for use with other integration methods.
    scRNA <- readRDS(file = "data/scRNA/ureter-scRNA.Rds")
    genes<-read.csv(file ="data/scRNA/genes.csv")[,1]

    # Source the functions file to make them available to our active session of R
    source("scripts/functions.R")

    ## ----- Load Visium samples ------

    U2.Seurat <- preProcessSeuratVisium("data/U2", normalization = "LogNormalize")

**Generate plots**

    #----- Perform anchor mapping----

    # This function will find "anchor" genes that are represented prominently in both the scRNA and SGE data. It will then use these to help approximate the distribution of single-cell identities in the spatial data.

    U2.Seurat <- anchorMapping(scRNA, U2.Seurat, feats = genes, query.dims=30, anchor.labels=levels(as.factor(scRNA$subclass)))

    # Plot the predictions made by the single-cell data

    scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(x) SpatialFeaturePlot(U2.Seurat, x))

    patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Seurat/Figure_1b.pdf", width =25, height = 25, units = "in", dpi = 300)

**OPTIONAL:** To generate plots with larger fonts use the following set
of code instead

    scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(x) SpatialFeaturePlot(U2.Seurat, x) +
                          theme(legend.key.size = unit(10, "mm"),
                                legend.text = element_text(size = 15),
                                legend.title = element_text(size = 20)))

    patchwork::wrap_plots(scplots, ncol=4) %T>% ggsave(filename = "figures/Seurat/Figure_1b.pdf", width =25, height = 25, units = "in", dpi = 300)

## Perform integration with Giotto

Again we will perform the integration only using `U2` Visium sample
provided via github. Please run the script `Giotto.R` or follow the code
and instructions below:

    # If running a new session (recommended), please run these two lines to load packages from your R library
    mylib="/home/sonas/R/mylib/" # replace this with your own R library path
    .libPaths(c(mylib, .libPaths()))

    library(ggplot2)
    library(Giotto)
    library(viridis)


    # This should be your directory/path where your sge-integration repository resides
    mydir='/home/sonas/star_protocol/'

    # Set this to your sge-integration folder
    setwd(paste0(mydir,"sge-integration/"))

    source("scripts/functions.R")



    scRNA <- readRDS("data/scRNA/ureter-scRNA.Rds")

    #----- Configure workspace with Giotto

    results_folder = paste0(mydir,'figures/Giotto')
    #mypython='/home/sonas/.local/share/r-miniconda/envs/giotto_env/bin/python3.6'
    instrs = createGiottoInstructions(save_dir = results_folder,
                                      save_plot = TRUE,
                                      show_plot = FALSE)

    #my_python_path = NULL


    ## Create new Giotto objects for each set of Visium data. We have performed the adjustment necessary, however new sets of Visium data will need to adjust the min_adj and max_adj values as necessary.
    # For more details, please see the Giotto vignette. 

    U2 <-  createGiottoVisiumObject(visium_dir = 'data/U2', expr_data = 'filter', 
                                         h5_visium_path = 'data/U2/filtered_feature_bc_matrix.h5',
                                         h5_tissue_positions_path = 'data/U2/spatial/tissue_positions_list.csv', 
                                         h5_image_png_path = 'data/U2/spatial/tissue_lowres_image.png',
                                    
                                         gene_column_index = 2, instructions = instrs, xmax_adj = 2000, ymin_adj = 1500, ymax_adj = 1600, xmin_adj = 1400)


    #These functions can help you determine whether your alignment is correct. See troubleshooting section of the STAR protocol for more details. 
    # spatPlot(gobject = U1, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
    #          save_param = list(save_name = 'U1_spatplot_image'))
    # 
    # spatPlot(gobject = U2, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
    #          save_param = list(save_name = 'U2_spatplot_image'))

    U2.Giotto<-preProcessGiotto(U2, "U2")

    # First you have to create the rank method signature matrix.
    sc_sign_matrix <- makeSignMatrixRank(
      sc_matrix = as.matrix(scRNA@assays$RNA@data),
      sc_cluster_ids = scRNA$subclass,
      ties_method = c("random"),
      gobject = NULL
    )

    #sc_sign_matrix<-as.matrix(sc_sign_matrix)

    U2.Giotto <- runSpatialEnrich(
      U2.Giotto,
      enrich_method ="rank",
      sign_matrix = sc_sign_matrix,
      expression_values = "normalized",
    )

    # Plot single cell enrichment patterns and simultaneously save figures in folder.
    scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = U2.Giotto, cell_color=unlist(c(U2.Giotto@spatial_enrichment$rank[,..i])), point_size = 2) + ggtitle(i) + scale_fill_distiller(palette = "Spectral"))

    patchwork::wrap_plots(scplots, ncol=4)%T>% ggsave(filename = "figures/Giotto/Figure_2b.pdf", width = 25, height = 20, units = "in", dpi = 300)

**OPTIONAL:** Use the following code for printing larger fonts for
figures

    scplots <- purrr::map(levels(as.factor(scRNA$subclass)), function(i) spatPlot(gobject = U2.Giotto, cell_color=unlist(c(U2.Giotto@spatial_enrichment$rank[,..i])), point_size = 2) +  
    theme(title = element_text(size=18),
              legend.text = element_text(size = 15),
              legend.title = element_text(size = 15),
              axis.title = element_text(size=15),
              axis.text = element_text(size=15)) + 
    ggtitle(i) + scale_fill_distiller(palette = "Spectral"))

    patchwork::wrap_plots(scplots, ncol=4) %T% ggsave(filename = "figures/Giotto/Figure_2b.pdf", width = 25, height = 20, units = "in", dpi = 300)

## Perform deconvolution with SPOTlight

    # Set this to your sge-integration folder
    mydir='/home/sonas/star_protocol/'

    setwd(paste0(mydir,"sge-integration/"))

    source("scripts/functions.R")

    library(Matrix)
    library(data.table)
    library(SPOTlight)
    library(igraph)
    library(RColorBrewer)
      

    scRNA <- readRDS("data/scRNA/ureter-scRNA.Rds")

    #----- Load data

    # Since SPOTlight uses the same format as Seurat, we can repurpose the methods to load
    U2.SL <- preProcessSeuratVisium("data/U2", normalization = "LogNormalize")

    # You may need to generate your own markers. This will take much longer if you do we give the option of loading some from a CSV file for the sample data.
    #markers <- read.csv("data/scRNA/markers.csv")
    Idents(scRNA) <- scRNA$subclass
    markers <- FindAllMarkers(scRNA, group.by = "subclass")

    colnames(markers)<-gsub('cluster','subclass',colnames(markers))

    mod_ls <- trainNMF(
        x = scRNA,
        y = U2.SL,
        groups = scRNA$subclass,
        mgs = markers,n_top = 10,
        weight_id = "avg_log2FC",
        group_id = "subclass",
        gene_id = "gene")

    res <- runDeconvolution(
        x = spe,
        mod = mod_ls[["mod"]],
        ref = mod_ls[["topic"]])

    U2.SL <- spotlightDeconvolve(vis = U2.SL, scrna = scRNA,markers =  markers,hvg=genes)

    cell.names <- names(which(colSums(U2.SL@meta.data[,-c(1:4, ncol(U2.SL@meta.data))]) > 0))


    ggsave(plot =  patchwork::wrap_plots(purrr::map(cell.names, function(x) SPOTlight::spatial_scatterpie(se_obj = U2.SL,
        cell_types_all = cell.names,
         cell_types_interest=x, 
         img_path = "data/U2/spatial/tissue_lowres_image.png",
         pie_scale = 0.4)), ncol = 4), filename="figures/SPOTlight/Figure_3b.pdf", 
        device="pdf", width = 25, height = 25, units = "in", dpi = 300)

**OPTIONAL:** Added another version of plotting function to print larger
fonts in the figures

    cell.names <- names(which(colSums(U1.SL@meta.data[,-c(1:4, ncol(U1.SL@meta.data))]) > 0))
    scplots<-purrr::map(cell.names, function(x) SPOTlight::spatial_scatterpie(se_obj = U1.SL,  cell_types_all = cell.names,  cell_types_interest=x,  
    img_path = "data/U1/spatial/tissue_lowres_image.png",
    pie_scale = 0.4)+
     theme(title = element_text(size=20),
           legend.text = element_text(size = 18),
           legend.title = element_text(size = 18),
           axis.title = element_text(size=18),
           axis.title.y = element_text(angle = 90),
           axis.text = element_text(size=18)))

    patchwork::wrap_plots(scplots, ncol=3) %T>% ggsave(filename="figures/SPOTlight/Figure3a.pdf", device="pdf", width = 25, height = 30, units = "in", dpi = 300)


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

    sessionInfo()

    ## R version 4.1.3 (2022-03-10)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: CentOS Linux 7 (Core)
    ## 
    ## Matrix products: default
    ## BLAS:   /opt/R-4.1.3/lib64/R/lib/libRblas.so
    ## LAPACK: /opt/R-4.1.3/lib64/R/lib/libRlapack.so
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_4.1.3  magrittr_2.0.3  fastmap_1.1.0   cli_3.4.1      
    ##  [5] tools_4.1.3     htmltools_0.5.3 rstudioapi_0.14 yaml_2.3.6     
    ##  [9] stringi_1.7.8   rmarkdown_2.17  knitr_1.40      stringr_1.4.1  
    ## [13] xfun_0.34       digest_0.6.30   rlang_1.0.6     evaluate_0.17
