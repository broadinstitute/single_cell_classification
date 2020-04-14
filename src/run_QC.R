source(here::here('src', 'load_data_helpers.R'))
source(here::here('src', 'QC_methods.R'))


run_all_QC <- function() {
  # run QC methods
  dat <- load_sc_exp(here::here('data'))
  seuObj <- create_seurat_obj(dat)
  num_CLs <- length(unique(seuObj@meta.data$singlet_ID))
  
  # set parameters based on the number of cell lines in the experiment
  n_pcs <- num_CLs*2
  param_range <- 10
  if(num_CLs > 25-param_range & num_CLs < 25+param_range) {
    clust_res <- 1
  } else if(num_CLs > 50-param_range & num_CLs < 50+param_range) {
    clust_res <- 2
  } else if(num_CLs > 100-param_range & num_CLs < 100+param_range) {
    clust_res <- 4
  } else {
    stop('Not implemented')
  }
  
  seuObj <- filter_low_quality_cells(seuObj)
  seuObj <- process_Seurat_obj(seuObj, n_pcs = n_pcs,  clust_res = clust_res)
  seuObj <- empty_droplet_clusters(seuObj)
  seuObj <- classify_doublets(seuObj)
  seuObj <- identify_low_confidence_cells(seuObj)
  
  return(seuObj)
}
