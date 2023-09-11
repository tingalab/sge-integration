## Setup 

It may be more intuitibe or efficient to download file using linux commands instead of R commonds. Use the following code do so:


`homedir='/path/to/your/home/dir/'
cd $homedir/sge-integration/

#download spatial data and scRNA supplemental files**
wget https://github.com/tingalab/sge-integration/releases/download/V1/protocol-data.tar.gz

#unzip
protocol-data.tar.gz 

#untar
tar -xvf protocol-data.tar 

#Remove compressed file
rm protocol-data.tar 

#Download scRNA data into the "data/scRNA/" folder
cd data/scRNA/
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE184nnn/GSE184111/suppl/GSE184111_scRNA_ureter10_normalized_counts.txt.gz

#Uncompress
gunzip GSE184111_scRNA_ureter10_normalized_counts.txt.gz

#Rename
mv GSE184111_scRNA_ureter10_normalized_counts.txt count-matrix.txt`

## Functions

Sometimes the function `Load10X_Spatial` generates error even when Seurat is loaded and all dependencies exist with proper data organization. To resolve this, simply add the following function code in the `functions.R` and source the fuctions script again or simply load the following function in the global environment.


Load10X_Spatial <- function(
    data.dir,
    filename = 'filtered_feature_bc_matrix.h5',
    assay = 'Spatial',
    slice = 'slice1',
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = NULL,
    ...
) {
  if (length(x = data.dir) > 1) {
    warning("'Load10X_Spatial' accepts only one 'data.dir'", immediate. = TRUE)
    data.dir <- data.dir[1]
  }
  data <- Read10X_h5(filename = file.path(data.dir, filename), ...)
  if (to.upper) {
    rownames(x = data) <- toupper(x = rownames(x = data))
  }
  object <- CreateSeuratObject(counts = data, assay = assay)
  if (is.null(x = image)) {
    image <- Read10X_Image(
      image.dir = file.path(data.dir, 'spatial'),
      filter.matrix = filter.matrix
    )
  } else {
    if (!inherits(x = image, what = "VisiumV1"))
      stop("Image must be an object of class 'VisiumV1'.")
  }
  image <- image[Cells(x = object)]
  DefaultAssay(object = image) <- assay
  object[[slice]] <- image
  return(object)
}


## Giotto related

If the suggested Giotto installation doesn't work for you, please install the alternate version using the following code:

devtools::install_github("drieslab/Giotto@suite")

It is imperative to use the updated the Giotto function code in `functions.R` script, if you switch to the above version:

preProcessGiotto2<-function(gobject, name){
  metadata = pDataDT(gobject)
  in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
  gobject <- subsetGiotto(gobject, cell_ids = in_tissue_barcodes)
  gobject <- normalizeGiotto(gobject, scalefactor = 10000, verbose = T)
  gobject <- addStatistics(gobject)
  gobject <- calculateHVF(gobject, save_param = list(save_name = paste0(name, '_HVFplot')),HVFname = "hvf",)
  gene_metadata = fDataDT(gobject)
  featgenes = gene_metadata[hvf == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID
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



This modified function is already added to the `fucntions.R` script. Simply update line 42 in the `Giotto.R` script to:

U2.Giotto<-preProcessGiotto2(U2, "U2")

