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

# this method run methods to classify each single cell based on reference profiles
# it assumes that the input data is stored in a folder data instead the single_cell_classification folder
# The required input data is:
# bulk_ref, bulk_alt: csv matrices of genes x samples for alt and ref allele counts for the reference profiles
# sc_ref, sc_alt: sparse matrices of single cells x genes of alt and ref allele counts for the single cell data, 
# column order must match row order from the reference profiles
# barcodes.tsv: tsv, single column file of just cell barcodes (with no column name)
run_SNP_classification <- function(run_stats = F, fast = F, call_doublets = T) {

  bulk_alt_p <- load_reference_files(here::here('data', 'bulk_ref.csv'), here::here('data', 'bulk_alt.csv'))

  sc_mats <- load_single_cell_files(here::here('data', 'sc_ref.mtx'), here::here('data', 'sc_alt.mtx'))

  cell_ids <- readr::read_tsv(here::here('data', 'barcodes.tsv'), col_names = F) %>% as.data.frame()

  # create output matrix
  num_cells <- nrow(sc_mats$sc_alt)

  SNP_class <- create_classification_matrix(num_cells, call_doublets, fast, run_stats)

  for(id in 1:num_cells) {
    cur_cell_res <- sc_SNP_match(sc_mats$sc_ref[id,,drop=FALSE], sc_mats$sc_alt[id,,drop=FALSE], bulk_alt_p, max_retry=5, doublet=call_doublets, ind=id, all_stats = run_stats, fast=fast)
    if(run_stats) {
      start_ind <- ((id-1) * ncol(bulk_alt_p)) + 1
      end_ind <- start_ind + (ncol(bulk_alt_p)-1)
      SNP_class[start_ind:end_ind,] <- cur_cell_res
    } else {
      SNP_class[id,] <- cur_cell_res
    }
  }

  colnames(SNP_class) <- colnames(cur_cell_res)
  SNP_class$barcode <- cell_ids[,1]
  
  return(SNP_class)
}

