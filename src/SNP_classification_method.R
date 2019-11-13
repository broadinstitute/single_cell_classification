
source(here::here('src', 'SNP_classification_helpers.R'))

sc_SNP_match <- function(ref_reads, alt_reads, bulk_alt_probs, max_retry = 5, doublet, ind, all_stats=F, fast=F) {
  
  #extract SNPs where we have at least some reads in the cell
  used_SNPs <- which(Matrix::colSums(ref_reads + alt_reads) > 0)
  
  if(length(used_SNPs) <= 1) {
    print("NO READS")
    # if no reads at any of the SNP sites return NA
    res <- return_na_singlet(ind, bulk_alt_probs, stats=all_stats, margin_stats=T)
    
  } else {
    # subset to SNP sites with reads
    ref_reads <- ref_reads[,used_SNPs]
    alt_reads <- alt_reads[,used_SNPs]
    
    sc_obs_mat <- cbind(ref_reads, alt_reads)
    
    tot_reads <- sum(ref_reads + alt_reads)
    num_SNPs <- length(used_SNPs)
    
    if(fast) {
      # doesn't produce margin stats, but runs faster and is used for out of pool classification
      res <- single_cell_classification(bulk_alt_probs[used_SNPs,], sc_obs_mat, num_coefs=1, max_retry=max_retry)
      
    } else {
      # produce margin stats which give confidence metrics of classification
      res <- singlet_classification_with_margin(bulk_alt_probs[used_SNPs,], sc_obs_mat, stats=all_stats, ind=ind)
    }
    res$tot_reads <- tot_reads
    res$num_SNPs <- num_SNPs
    
    #find best doublet pair
    if(doublet==TRUE) {
      top_singlet <- res$singlet_ID %>% as.character()
      doublet_class <- single_cell_classification(bulk_alt_probs[used_SNPs,], sc_obs_mat, num_coefs=2, max_retry=max_retry, top_singlet =top_singlet, singlet_dev = res$singlet_dev)
      res$doublet_ratio <- doublet_class$doublet_ratio
      if(is.na(doublet_class$doublet_dev_ratio) | is.na(res$singlet_dev)) {
        res$doublet_dev_imp <- NA
      } else {
        res$doublet_dev_imp <- doublet_class$doublet_dev_ratio - res$singlet_dev
      }
      res$doublet_CL1 <- doublet_class$doublet_CL1
      res$doublet_CL2 <- doublet_class$doublet_CL2
    }
  }
  
  res$singlet_ID <- as.character(res$singlet_ID)
  return(res)
}

