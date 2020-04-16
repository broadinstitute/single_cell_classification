options(repos = c("https://cran.cnr.berkeley.edu"))

cran_packages <- c('here', 'magrittr', 'reshape2', 'plyr', 'dplyr', 'ggplot2',
                   'Matrix', 'readr', 'tibble', 'data.table', 'glmnet', 'Seurat',
                   'vcfR', 'mclust')
# package versions: here (v 0.1), magrittr (v 1.5), reshape2 (v 1.4.3), plyr (v1.8.5), dplyr (v 0.8.3),
# ggplot2 (v 3.2.1), Matrix (v 1.2.18), readr (v 1.3.1), tibble (v 2.1.3), data.table (v 1.12.8),
# glmnet (v 3.0.2), Seurat (v.3.1.2), vcfR (v 1.9.0), mclust (v 5.4.5)

new_cran_packages <- cran_packages[!(cran_packages %in% installed.packages()[,"Package"])]
if(length(new_cran_packages)) install.packages(new_cran_packages)

bioconductor_packages <- c('DropletUtils', 'edgeR')
# package versions: DropletUtils (v 1.6.1), edgeR (v 3.28.0)
new_bioconductor_packages <- bioconductor_packages[!(bioconductor_packages %in% installed.packages()[,"Package"])]
if(length(new_bioconductor_packages)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  }
  BiocManager::install(new_bioconductor_packages)
}

