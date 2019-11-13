library(mclust)
library(ggplot)

# begin adding cell quality classifications
# classify cells with too few SNPs or abnormal mitochondrial read percentages as low quality
filter_low_quality_cells <- function(seuObj, min_SNPs = 50, max_mito_perc = 25, min_mito_perc=1) {
  
  #too few SNPs
  seuObj@meta.data <- seuObj@meta.data %>%
    dplyr::mutate(cell_quality = ifelse(num_SNPs < min_SNPs, 'low_quality', 'normal'))
  
  #abnormal mitochondrial read %
  seuObj@meta.data <- seuObj@meta.data %>%
    dplyr::mutate(cell_quality = ifelse(percent.mito > max_mito_perc | percent.mito < min_mito_perc, 'low_quality', cell_quality))
  
  return(seuObj)
}

# run Seurat methods, excluding low quality cells
process_Seurat_obj <- function(seuObj, n_pcs = 50, cluster_k =10, clust_res = 1, plot = TRUE) {
  
  all_meta <- seuObj@meta.data
  rownames(all_meta) <- all_meta$barcode
  
  rownames(seuObj@meta.data) <- seuObj@meta.data$barcode
  
  seuObj <- Suerat::subset(seuObj, cells = dplyr::filter(seuObj@meta.data, cell_quality == 'normal')$barcode)
  # run Seurat methods to get gene expression clusters
  seuObj <- Seurat::NormalizeData(object = seuObj, 
                                  normalization.method = "LogNormalize", 
                                  scale.factor = 1e5, verbose = FALSE)
  
  seuObj <- Seurat::ScaleData(object = seuObj, vars.to.regress = c(), verbose = FALSE)
  
  seuObj <- Seurat::FindVariableFeatures(object = seuObj,
                                         nfeatures = 5000,
                                         do.plot = FALSE,
                                         selection.method = 'vst', verbose = FALSE)
  
  seuObj <- Seurat::RunPCA(object = seuObj,
                           features = Seurat::VariableFeatures(seuObj),
                           seed.use = 1,
                           npcs = n_pcs,
                           verbose = FALSE)
  
  seuObj <- Seurat::RunTSNE(object = seuObj,
                            dims = 1:n_pcs,
                            check_duplicates = FALSE,
                            seed.use = 1,
                            perplexity = 30,
                            verbose = FALSE)  
  
  seuObj <- Seurat::FindNeighbors(seuObj, reduction = 'pca',
                                  dims = 1:n_pcs,
                                  k.param = cluster_k, 
                                  force.recalc = TRUE,
                                  verbose = FALSE)
  seuObj <- Seurat::FindClusters(seuObj, resolution = clust_res, random.seed = 1, verbose = FALSE)
  
  if(plot) {
    gg <- Seurat::DimPlot(seuObj, group.by = 'seurat_clusters') + Seurat::NoLegend()
    print(gg)
  }
  
  all_meta <- dplyr::left_join(all_meta, seuObj@meta.data[,c('barcode', 'RNA_snn_res.1','seurat_clusters')])
  seuObj@meta.data <- all_meta
  return(seuObj)
  
}

# identify low quality 'cells', which may be from empty droplets
empty_droplet_clusters <- function(seuObj, dev_thresh = 0.3, plot = TRUE) {
  df <- Seurat::Embeddings(seuObj, reduction = 'tsne') %>% as.data.frame() %>%
    tibble::rownames_to_column(var = 'barcode') %>%
    dplyr::left_join(seuObj@meta.data %>% 
                       dplyr::mutate(max_dev = singlet_dev + doublet_dev_imp), by='barcode')
  
  #compute stats for each cluster
  clust_stats <- plyr::ddply(df, .(seurat_clusters), function(ss) {
    data_frame(n_cells = nrow(ss),
               med_doub = median(ss$doublet_dev_imp, na.rm=T),
               med_dev = median(ss$max_dev, na.rm=T))
  }) %>% dplyr::mutate(
    clust_type = ifelse(med_dev <= dev_thresh, 'empty_droplet', 'normal')
  )
  
  empty_droplets <- df$barcode[which((clust_stats$clust_type %>% rlang::set_names(clust_stats$seurat_clusters) %>%
                                        .[df$seurat_clusters]) == 'empty_droplet')]
  print(sprintf('Identified %d empty droplets', length(empty_droplets)))
  
  rownames(seuObj@meta.data) <- seuObj@meta.data$barcode
  seuObj@meta.data[empty_droplets,'cell_quality'] <-  'empty_droplet'
  
  if(plot) {
    gg <- Seurat::DimPlot(seuObj, group.by = 'cell_quality')
    print(gg)
  }
  return(seuObj)
}

# Identify doublets
classify_doublets <- function(seuObj, doub_prob_thresh=0.5, doublet_dev_offset = 0.001, plot = TRUE) {
  # use singlet deviance, doublet deviance improvement, and fraction of genes detected in a given cell
  # to identify cells that are doublets
  # doublet identification
  
  doublet_df <- Seurat::FetchData(seuObj, c('singlet_dev', 'doublet_dev_imp', 'cell_quality')) %>% 
    tibble::rownames_to_column(var = 'barcode') %>% 
    dplyr::filter(cell_quality %in% c('normal', 'doublet')) %>% 
    dplyr::mutate(log_doublet_imp = log10(doublet_dev_imp + doublet_dev_offset)) %>% 
    dplyr::select(-doublet_dev_imp) 
  
  
  # handle cells where the doublet model failed by treating them as singlets.
  fail_doublet_model <- doublet_df %>% 
    dplyr::filter(is.na(log_doublet_imp)) %>% 
    .[['barcode']]
  print(sprintf('%d cells with failed doublet models treated as singlets', length(fail_doublet_model)))
  
  doublet_df %<>% dplyr::filter(!(barcode %in% fail_doublet_model))
  
  # sfit GMM
  gmm_output <- mclust::Mclust(as.matrix(doublet_df[, c('log_doublet_imp', 'singlet_dev')]), G = 2, verbose = FALSE, prior = mclust::priorControl(shrinkage = 0))
  
  doublet_cluster <- which.max(gmm_output$parameters$mean['log_doublet_imp', ])
  singlet_cluster <- which.min(gmm_output$parameters$mean['log_doublet_imp', ])
  doublet_df$doublet_prob <- gmm_output$z[,doublet_cluster]
  doublet_cells <- doublet_df %>% 
    dplyr::filter(doublet_prob >= doub_prob_thresh) %>% 
    .[['barcode']]
  
  
  seuObj@meta.data[doublet_cells, 'cell_quality'] <- 'doublet'
  
  if(plot) {
    g1 <- plot(gmm_output, what = 'classification')
    
    g2 <- ggplot2::ggplot(doublet_df, ggplot2::aes(singlet_dev, log_doublet_imp, color=doublet_prob)) + 
      ggplot2::geom_point() + 
      ggplot2::theme(legend.position = "right", legend.direction = "vertical", legend.text = ggplot2::element_text(size=10)) +
      ggplot2::xlab("singlet deviance") + ggplot2::ylab("log doublet deviance improvement")
    
    
    g3 <- ggplot2::ggplot(doublet_df %>% dplyr::mutate(is_doublet = barcode %in% doublet_cells), ggplot2::aes(singlet_dev, log_doublet_imp, color= is_doublet)) +
      ggplot2::geom_point() +
      ggplot2::theme(legend.position = "right", legend.direction = "vertical", legend.text =ggplot2:: element_text(size=10)) +
      ggplot2::xlab("singlet deviance") + ggplot2::ylab("log doublet deviance improvement") +
      ggplot2::guides(color = ggplot2::guide_legend(title = 'doublet'))
    
    print(g1)
    print(g2)
    print(g3)
  }
  
  return(seuObj)
}

identify_low_confidence_cells <- function(seuObj, min_z_margin = 2, plot = TRUE) {
  df <- seuObj@meta.data
  
  seuObj@meta.data %<>% 
    dplyr::mutate(cell_quality = ifelse(seuObj@meta.data$singlet_z_margin < min_z_margin & seuObj@meta.data$cell_quality == 'normal', 'low_confidence', cell_quality))
  
  if(plot) {
    cell_class_confidence <- ggplot2::ggplot(df %>% filter(cell_quality == 'normal'),
                                             ggplot2::aes(singlet_z_margin)) + 
      ggplot2::geom_density() + 
      ggplot2::geom_vline(xintercept = min_z_margin, linetype = 'dashed')
    
    print(cell_class_confidence)
    gg <- ggplot2::ggplot(seuObj@meta.data, ggplot2::aes(singlet_dev, doublet_dev_imp, fill=cell_quality)) + 
      ggplot2::geom_point(pch = 21, size = 1.5, color = 'white', stroke = 0.1) + ggplot2::ggtitle('cell quality classifications')  
    print(gg)
  }
  
  return(seuObj)
}

# run QC methods
run_QC <- function(data_folder) {
  dat <- load_sc_exp(data_folder)
  seuObj <- create_seurat_obj(dat)
  num_CLs <- length(unique(seuObj@meta.data$singlet_ID))
  n_pcs <- num_CLs*2
  if(num_CLs < 40) {
    clust_res <- 1
  } else if(num_CLs>40 & num_CLs <60) {
    clust_res <- 2
  } else if(num_CLs>90 & num_CLs <110) {
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
