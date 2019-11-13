source(here::here('src', 'QC_methods.R'))



# run QC methods
dat <- load_sc_exp(here::here('data'))
seuObj <- create_seurat_obj(dat)
low_quality_cells <- identify_low_quality_cells(seuObj)
doublets <- classify_doublets(seuObj, low_quality_cells)
seuObj <- add_cell_quality_classifications(seuObj, low_quality_cells, doublets)

