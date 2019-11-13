# return NA for cells with no reads at the given SNP sites
return_na_singlet <- function(cur_ind, bulk_alt_probs, stats=F, margin_stats = F) {
  if(margin_stats) {
    df <- data.frame(singlet_ID = NA,
                     singlet_dev = NA,
                     singlet_dev_z = NA,
                     singlet_margin = NA,
                     singlet_z_margin = NA,
                     doublet_z_margin = NA)
  } else {
    df <- data.frame(singlet_ID = NA,
                     singlet_dev = NA,
                     deviance = NA)
  }
  
  if(stats) {
    df <- as.data.frame(matrix(data=0, nrow=ncol(bulk_alt_probs), ncol=5))
    colnames(df) <- c("singlet_ID", "coef", "dev_ratio", "dev_ratio_z", "singlet_deviance")
    df$CCLE_ID <- colnames(bulk_alt_probs)
    df$barcode <- rep(cur_ind, ncol(bulk_alt_probs))
  }
  
  return(df)
}

# classification model used for doublet classifications (can also be used for faster singlet classifications w/out stats)
single_cell_classification <- function(bulk_alt_probs, sc_obs_mat, num_coefs, max_retry, top_singlet=NA, singlet_dev = 0) {
  if(colSums(sc_obs_mat)[1] == 0 | colSums(sc_obs_mat)[2] == 0) {
    res <- list(doublet_dev_ratio = NA, doublet_CL1 = NA, doublet_CL2=NA)
    return(res)
  }
  mod <- glmnet::glmnet(x = bulk_alt_probs, 
                        y = sc_obs_mat, 
                        family = 'binomial',
                        alpha = 1,
                        dfmax = num_coefs,
                        lower.limits = 0,
                        nlambda = 100)
  
  mod_devs <- deviance(mod)
  best_model <- which(mod$df == num_coefs) %>% rev() %>% head(1)
  
  retry_count <- 0
  while(length(best_model) == 0 & retry_count <= max_retry &
        !is.infinite(mod$lambda %>% rev() %>% head(1)) &
        length(mod$lambda[which(mod$df < 2) %>% rev() %>% head(1)]) != 0) {
    
    too_small <- mod$lambda %>% rev() %>% head(1)
    too_big <- mod$lambda[which(mod$df < num_coefs) %>% rev() %>% head(1)]
    lambdas <- seq(from = too_big, to = too_small, length.out = 100)
    mod <- glmnet::glmnet(x = bulk_alt_probs, 
                          y = sc_obs_mat, 
                          family = 'binomial',
                          alpha = 1,
                          dfmax = num_coefs,
                          lower.limits = 0,
                          lambda = lambdas)
    mod_devs <- deviance(mod)
    best_model <- which(mod$df == num_coefs) %>% rev() %>% head(1)
    retry_count <- retry_count + 1
  }
  
  if (length(best_model) == 0 | length(unique(mod$lambda))==1) {
    if(num_coefs  == 1) {
      res <- return_na_singlet(ind, bulk_alt_probs)
    } else {
      # if model doesn't produce 2 non-zero cell lines test individually w/ cell lines
      res <- run_doublet_model_per_CL(top_singlet, bulk_alt_probs, sc_obs_mat)
    }
  } else {
    coefs <- coef(mod, s = mod$lambda[best_model])[2:(sum(ncol(bulk_alt_probs))+1),1] %>% 
      sort(decreasing = TRUE)
    
    CL1 <- names(coefs)[1]
    if(num_coefs==1) {
      bulk_profile <- bulk_alt_probs[,CL1]
      
      #fit singlet model without regularization to get stats
      glm_input <- cbind.data.frame(sc_obs_mat, bulk_profile)
      colnames(glm_input) <- c('sc_ref', 'sc_alt', 'bulk_CL1')
      glm_mod <- glm(cbind(sc_ref, sc_alt) ~ ., data = glm_input, family = binomial("logit")) %>% summary()
      
      null.deviance <- glm_mod$null.deviance
      deviance <- glm_mod$deviance
      dev_ratio <- 1-(deviance/null.deviance)
      res <- data.frame(
        deviance = deviance,
        singlet_dev = dev_ratio,
        singlet_ID = CL1
      )
      
    } else {
      CL2 <- names(coefs)[2]
      cur_CLs <- c(CL1, CL2)
      
      #fit doublet model without regularization to get stats
      glm_input <- cbind.data.frame(sc_obs_mat, bulk_alt_probs[,cur_CLs])
      colnames(glm_input) <- c('sc_ref', 'sc_alt', 'bulk_CL1', 'bulk_CL2')
      glm_mod <- glm(cbind(sc_ref, sc_alt) ~ ., data = glm_input, family = binomial("logit")) %>% summary()
      
      null.deviance <- glm_mod$null.deviance
      deviance <- glm_mod$deviance
      doublet_dev_ratio <- 1-(deviance/null.deviance)
      
      if(!is.na(doublet_dev_ratio) & !is.na(singlet_dev) & doublet_dev_ratio < singlet_dev) {
        res <- run_doublet_model_per_CL(top_singlet, bulk_alt_probs, sc_obs_mat)
        doublet_dev_ratio <- res$doublet_dev_ratio
        CL1 <- res$doublet_CL1
        CL2 <- res$doublet_CL2
      }
      res <- list(doublet_dev_ratio = doublet_dev_ratio, doublet_CL1 = CL1, doublet_CL2=CL2)
    }
  }
  
  return(res)
}

# take the top cell line from singlet model and test with all other cell lines
# to find most likely second cell line if a doublet
run_doublet_model_per_CL <- function(top_singlet, bulk_alt_probs, sc_obs_mat) {
  top_prof <- bulk_alt_probs[,top_singlet, drop=F]
  other_bulk_probs <- bulk_alt_probs[,setdiff(colnames(bulk_alt_probs), top_singlet)]
  ind_bulk_matches <- adply(other_bulk_probs, 2, function(cur_bulk_profile) {
    #fit doublet model without regularization to get stats
    glm_input <- cbind.data.frame(sc_obs_mat, top_prof, cur_bulk_profile)
    colnames(glm_input) <- c('sc_ref', 'sc_alt', 'bulk_CL1', 'bulk_CL2')
    glm_mod <- glm(cbind(sc_ref, sc_alt) ~ ., data = glm_input, family = binomial("logit")) %>% summary()
    null.deviance <- glm_mod$null.deviance
    deviance <- glm_mod$deviance
    doublet_dev_ratio <- 1-(deviance/null.deviance)
    
    data.frame(
      dev_ratio = doublet_dev_ratio
    )
  }, .id = 'CCLE_ID')
  
  ind_bulk_matches %<>% 
    dplyr::arrange(dplyr::desc(dev_ratio))
  
  return(list(doublet_dev_ratio=ind_bulk_matches$dev_ratio[1], 
              doublet_CL1=top_singlet, doublet_CL2=as.character(ind_bulk_matches$CCLE_ID[1])))
}

# identify best matching reference cell line for each cell
singlet_classification_with_margin <- function(bulk_alt_probs, sc_obs_mat, stats, tot_reads, num_SNPs, ind){
  
  # calculate individually for each cell line in order to get
  # stat about the deviance margin
  
  # calc glm for each bulk cell line
  ind_bulk_matches <- plyr::adply(bulk_alt_probs, 2, function(bulk_profile) {
    glm_input <- cbind.data.frame(sc_obs_mat, bulk_profile)
    colnames(glm_input) <- c('sc_ref', 'sc_alt', 'bulk_frac')
    glm_mod <- glm(cbind(sc_ref, sc_alt) ~ bulk_frac, data = glm_input, family = binomial("logit")) %>% summary()
    
    data.frame(
      null.deviance = glm_mod$null.deviance,
      deviance = glm_mod$deviance,
      num_SNPs = length(bulk_profile)
    )
  }, .id = 'singlet_ID')
  
  # select cell line with lowest deviance
  ind_bulk_matches <- ind_bulk_matches %>% dplyr::arrange(deviance) 
  ind_bulk_matches$singlet_dev <- 1-(ind_bulk_matches$deviance/ind_bulk_matches$null.deviance)
  cur_sd <- max(sd(ind_bulk_matches$singlet_dev[-1], na.rm=T), 0.0001)
  # remove top cell line when calculated sd to account for differences in sd behavior due to different pool sizes
  ind_bulk_matches$dev_ratio_z <- (ind_bulk_matches$singlet_dev - mean(ind_bulk_matches$singlet_dev, na.rm=T)) / 
    cur_sd
  
  classification_res <- data.frame(singlet_ID = as.character(ind_bulk_matches$singlet_ID[1]),
                                   num_SNPs = ind_bulk_matches$num_SNPs[1],
                                   singlet_dev = ind_bulk_matches$singlet_dev[1],
                                   singlet_dev_z = ind_bulk_matches$dev_ratio_z[1],
                                   singlet_margin =  ind_bulk_matches$singlet_dev[1] - ind_bulk_matches$singlet_dev[2],
                                   singlet_z_margin = ind_bulk_matches$dev_ratio_z[1] - ind_bulk_matches$dev_ratio_z[2],
                                   doublet_z_margin = ind_bulk_matches$dev_ratio_z[2] - ind_bulk_matches$dev_ratio_z[3])
  
  if(stats) {
    return(ind_bulk_matches)
  } else {
    return(classification_res)
  }
}
