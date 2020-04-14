options(repos = c("https://cran.cnr.berkeley.edu"))

cran_packages <- c('here', 'magrittr', 'reshape2', 'plyr', 'dplyr', 'ggplot2',
                   'Matrix', 'readr', 'tibble', 'data.table', 'glmnet', 'Seurat',
                   'vcfR', 'mclust')
new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

bioconductor_packages <- c('DropletUtils', 'edgeR')
new_bioconductor_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new_bioconductor_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  }
  BiocManager::install(new_bioconductor_packages)
}

