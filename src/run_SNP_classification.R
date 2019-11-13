# run single cell classification
# identify best fitting reference sample for each cell using SNPs

library(plyr)
library(dplyr)
library(magrittr)
library(reshape2)
library(tibble)
library(readr)
library(glmnet)
library(here)

source(here::here('src', 'load_data_helpers.R'))
source(here::here('src', 'SNP_classification_helpers.R'))
source(here::here('src', 'SNP_classification_method.R'))

bulk_alt_p <- load_reference_files(here::here('data', 'bulk_ref.csv'), here::here('data', 'bulk_alt.csv'))

sc_mats <- load_single_cell_files(here::here('data', 'sc_ref.mtx'), here::here('data', 'sc_alt.mtx'))

cell_ids <- readr::read_tsv(here::here('data', 'barcodes.tsv'), col_names = F) %>% as.data.frame()

SNP_class <- plyr::ldply(seq(nrow(cell_ids)), function(id) {
  sc_SNP_match(sc_mats$sc_ref[,id], sc_mats$sc_alt[,id], bulk_alt_p, max_retry=5, doublet=TRUE, ind=id)
}) 

SNP_class$barcode <- cell_ids[,1] 