# run single cell classification, 
# downsampling the single cell alt and ref allele count matrices
# to evaluate our classification

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

#params
run_stats <- F
fast <- T
call_doublets <- F

bulk_alt_p <- load_reference_files(here::here('data', 'bulk_ref.csv'), here::here('data', 'bulk_alt.csv'))

sc_mats <- load_single_cell_files(here::here('data', 'sc_ref.csv'), here::here('data', 'sc_alt.csv'))

cell_ids <- readr::read_tsv(here::here('data', 'barcodes.tsv'), col_names = F)

# create output matrix
num_cells <- nrow(sc_mats$sc_alt)

sampling_rates = c(0.9, 0.8, 0.7, 0.6,0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.04, 0.03, 0.02, 0.01)
for(ds in sampling_rates) {
  downsampled_mats <- downsample_single_cell(sc_mats$sc_alt, sc_mats$sc_ref, ds) 
  SNP_class <- create_classification_matrix(num_cells, call_doublets, fast, run_stats)
  
  for(id in 1:num_cells) {
    cur_cell_res <- sc_SNP_match(downsampled_mats$sc_ref[id,,drop=FALSE], downsampled_mats$sc_alt[id,,drop=FALSE],
                                 bulk_alt_p, max_retry=5, doublet=call_doublets, ind=id, all_stats = run_stats, fast=fast)
    SNP_class[id,] <- cur_cell_res
  }
  
  colnames(SNP_class) <- colnames(cur_cell_res)
  SNP_class$barcode <- cell_ids
  readr::write_csv(SNP_class, here::here('data', paste0('classifications_',ds,'.csv')))
}
