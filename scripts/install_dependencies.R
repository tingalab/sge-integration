
#Install particular version of CRAN package 

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

# Now install other packages (Giotto)

# For detailed installation instructions and troubleshooting please see https://giottosuite.readthedocs.io/en/latest/gettingstarted.html

#remotes::install_github("drieslab/Giotto@suite")
remotes::install_github("RubD/Giotto@v1.1.0") 


# Run this once to setup Giotto environment

# installGiottoEnvironment(packages_to_install = c("pandas==1.1.5", "networkx==2.6.3", "python-igraph==0.9.6", "leidenalg==0.8.7","python-louvain==0.15", "scikit-learn==0.24.2"), force_environment = TRUE)
#            
 installGiottoEnvironment(packages_to_install = c("pandas", 
                                                 "networkx", 
                                                  "python-igraph", 
                                                  "leidenalg",
                                                  "python-louvain",
                                                  "scikit-learn"), 
                          force_environment = TRUE)

