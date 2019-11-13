# reference gene expression should be read level data, with rows as genes and columns as samples
# CCLE_RNAseq_reads taken from: https://depmap.org/portal/download/
# single cell input is taken from running load_sc_exp
transform_GE_data <- function(data, reference_GE) {
  sc_GE <- scale(as.matrix(data$counts), center = FALSE, scale = colSums(as.matrix(data$counts))) * 1e5
  sc_GE <- log2(1 + sc_GE) # log transform
  
  rownames(reference_GE) <- str_match(rownames(reference_GE), '\\((ENSG.+)\\)')[,2]
  #colnames(reference_GE) <- arxspan.to.ccle(colnames(reference_GE))
  reference_GE <- scale(reference_GE, center = FALSE, scale = colSums(reference_GE, na.rm=T)) * 1e6
  reference_GE <- log2(reference_GE + 1) # log transform
  
  # mean center expression profiles
  rel_sc_GE <- scale(t(sc_GE), center = TRUE, scale = FALSE) %>% t()
  rel_GE <- scale(t(reference_GE), center = TRUE, scale = FALSE) %>% t()
  
  return(list(rel_sc_GE, rel_GE))
  
}

# identify subset of genes that are found in a larger fraction of cells or are highly expressed
get_used_genes <- function(data, min_avg_GE=1, ext_perc=0.98, min_frac_expressing = 0.5, min_ext_exp = 3) {
  sc_GE <- scale(as.matrix(data$counts), center = FALSE, scale = colSums(as.matrix(data$counts))) * 1e5
  sc_GE <- log2(1 + sc_GE)
  
  gene_df <- data$gene_info
  
  gene_df$num_det_CLs = rowSums(data$counts > 0)
  gene_df$extreme_expression = matrixStats::rowQuantiles(as.matrix(sc_GE), probs = ext_perc)
  gene_df$avg_expression <- rowMeans(sc_GE)
  
  if (is.null(min_ext_exp)) {
    min_ext_exp <- Inf
  }
  # select genes that are expressed in at least a given fraction of cells or are highly expressed
  if (!is.null(min_frac_expressing)) {
    min_det_CLs <- ncol(sc_GE) * min_frac_expressing
    gene_df$included <- gene_df$num_det_CLs > min_det_CLs | gene_df$extreme_expression > min_ext_exp
  } else {
    gene_df %<>% mutate(included = avg_expression >= min_avg_GE)
  }
  cat(sprintf('Using %d/%d genes\n', sum(gene_df$included), nrow(gene_df)))
  return(gene_df)
}

# calculate single cell gene expression profiles that are smoothed across neighboring cells
calc_smoothed_GE <- function(rel_sc_GE, gene_df, initial_dims = 50, perplexity = 30, dims=3, theta = 0.2) {
  # PCA
  pca <- prcomp(rel_sc_GE[gene_df$included,])
  
   # t-SNE (using pca for initialization to create stable t-SNE)
   tsne_res <- Rtsne(t(as.matrix(rel_sc_GE[gene_df$included,])),
                    initial_dims = initial_dims,
                    perplexity = perplexity,
                    dims=dims, theta = theta, Y_init=pca$rotation[,1:dims])
  
  svm_band <- kernlab::sigest(tsne_res$Y)[3]
  gk <- KRLS::gausskernel(X=tsne_res$Y, sigma=svm_band)
  
  # computing smoothed sc ge in parallel
  weighted_mean <- function(y, rs_ge) {
    return(apply(rs_ge,1, function(x) weighted.mean(x,y)))
  }
  reduced_sc_ge <- rel_sc_GE[gene_df$included,]
  cl <- parallel::makeCluster(mc <- getOption("cl.cores", 4))
  parallel::clusterExport(cl=cl, varlist=c("weighted_mean","gk", "reduced_sc_ge"))
  smoothed_sc_ge <- parallel::parApply(cl,gk,2, function(y) weighted_mean(y, reduced_sc_ge))
  parallel::stopCluster(cl) 
  
  return(smoothed_sc_ge)
  
}

# classify single cell as given reference samples using the smoothed expression profiles
call_GE_regession_classifier <- function(rel_sc_GE, gene_df, smoothed_sc_ge, rel_GE, target_pool=NULL) {
  if(is.null(target_pool)) {
    target_pool <- colnames(ref_GE)
  }
  GE_res <- plyr::ldply(seq(1,ncol(rel_sc_GE[gene_df$included,])), function(id) {
    print(id)
    sc_exp_match(rel_GE[gene_df$included,target_pool], rel_sc_GE[gene_df$included,id], smoothed_sc_ge[gene_df$included,id], ind=id)
  }) 
  
  return(GE_res)
  
}

# identify the best matching singlet reference profile for each cell,
# and fit a doublet model
sc_exp_match <- function(cur_GE, cur_sc_GE, smoothed_sc_GE, max_retry=5, ind) {
  ind_bulk_matches <- plyr::adply(cur_GE, 2, function(bulk_profile) {
    cur_cor <- cor(bulk_profile, smoothed_sc_GE, use="pairwise")
    reg <- lm(cur_sc_GE ~ bulk_profile)
    data.frame(
      coef = cur_cor,
      dev_ratio = summary(reg)$r.squared
    )
  }, .id = 'CCLE_ID')
  
  ind_bulk_matches$dev_ratio_z <- (ind_bulk_matches$dev_ratio - mean(ind_bulk_matches$dev_ratio, na.rm=T)) / 
    sd(ind_bulk_matches$dev_ratio, na.rm=T)
  
  ind_bulk_matches %<>%
    dplyr::arrange(dplyr::desc(coef))
  best_singlet <- ind_bulk_matches %>% 
    head(1)
  best_singlet_z_margin <- ind_bulk_matches[1, 'dev_ratio_z'] - ind_bulk_matches[2, 'dev_ratio_z']
  #find best doublet pair
  mod <- glmnet::glmnet(x = cur_GE, 
                y = cur_sc_GE, 
                alpha = 1,
                dfmax = 2,
                lower.limits = 0,
                nlambda = 100,
                family='gaussian') #?
  mod_devs <- deviance(mod)
  best_doublet_model <- which(mod$df == 2) %>% rev() %>% head(1)
  
  
  retry_count <- 0
  while(length(best_doublet_model) == 0 & retry_count <= max_retry & !is.infinite(mod$lambda)) {
    too_small <- mod$lambda %>% rev() %>% head(1)
    too_big <- mod$lambda[which(mod$df < 2) %>% rev() %>% head(1)]
    lambdas <- seq(from = too_big, to = too_small, length.out = 100)
    mod <- glmnet::glmnet(x = cur_GE, 
                  y = cur_sc_GE, 
                  alpha = 1,
                  dfmax = 2,
                  lower.limits = 0,
                  lambda = lambdas)
    mod_devs <- deviance(mod)
    best_doublet_model <- which(mod$df == 2) %>% rev() %>% head(1)
    retry_count <- retry_count + 1
  }
  if (length(best_doublet_model) == 0) {
    doublet_dev_ratio <- NA
    doublet_dev_imp <- NA
    doublet_ratio <- NA
    doublet_CL1 <- NA
    doublet_CL2 <- NA
  } else {
    doublet_coefs <- coef(mod, s = mod$lambda[best_doublet_model])[2:(sum(ncol(cur_GE))+1),1] %>% 
      sort(decreasing = TRUE)
    doublet_CL1 <- names(doublet_coefs)[1]
    doublet_CL2 <- names(doublet_coefs)[2]
    
    #fit doublet model without regularization to get stats
    dmod <- glmnet::glmnet(x = cur_GE[,c(as.character(doublet_CL1), as.character(doublet_CL2))], 
                   y = cur_sc_GE, 
                   alpha = 0,
                   lower.limits = 0,
                   lambda = c(0),
                   family='gaussian')
    doublet_coefs <- coef(dmod)[2:3, 1] %>% sort()
    doublet_ratio <- doublet_coefs[1]/(sum(doublet_coefs))
    doublet_dev_ratio <- dmod$dev.ratio
    doublet_dev_imp <- doublet_dev_ratio - best_singlet$dev_ratio
    
  }
  df <-  data.frame(singlet_ID = best_singlet$CCLE_ID,
                    singlet_dev = best_singlet$dev_ratio,
                    singlet_dev_z = best_singlet$dev_ratio_z,
                    singlet_z_margin = best_singlet_z_margin,
                    doublet_ratio = doublet_ratio,
                    doublet_dev_imp = doublet_dev_imp,
                    doublet_CL1 = doublet_CL1,
                    doublet_CL2 = doublet_CL2)
  
  
  
  if(!is.na(doublet_CL2)) {
    singlet_mod <- lm(cur_sc_GE ~ cur_GE[,best_singlet$CCLE_ID])
    doublet_mod <- lm(cur_sc_GE ~ cur_GE[,c(doublet_CL1, doublet_CL2)])
    anov <- anova(singlet_mod, doublet_mod, test="Chisq")
    df$anova <- anov$`Pr(>Chi)`[2]
  } else {
    df$anova <- NA
  }
  return(df) 
}

