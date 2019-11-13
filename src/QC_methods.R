library(mclust)
library(ggplot2)

# take in counts matrix, gene info, and cell classifications for a given experiment 
# (input is produced from load_sc_exp)
# and create seurat object
create_seurat_obj <- function(dat) {
  dat <- dat %>%
    convert_to_hugo_symbols() #convert genes to hugo symbols (and only keep unique hugo symbol genes)
  
  rownames(dat$cell_info) <- dat$cell_info$barcode
  
  #make into Seurat object
  seuObj <- Seurat::CreateSeuratObject(dat$counts, 
                               min.cells = 0,
                               min.features = 0,
                               meta.data = dat$cell_info)
  
  
  #make singlet classifications the cell identifiers
  seuObj <- Seurat::SetIdent(seuObj, value = seuObj@meta.data$singlet_ID)
  
  #store gene info here
  seuObj@misc <- dat$gene_info
  
  #add mitochondrial gene fraction
  mito.genes <- grep(pattern = "^MT-",
                     x = rownames(seuObj@assays$RNA@counts),
                     value = TRUE) 
  percent.mito <- Matrix::colSums(seuObj@assays$RNA@counts[mito.genes, ])/Matrix::colSums(seuObj@assays$RNA@counts)
  seuObj@meta.data$percent.mito <- percent.mito
  
  #add cellular detection rate (fraction of genes detected in a given cell)
  seuObj[['cell_det_rate']] <- seuObj$nFeature_RNA/nrow(seuObj@assays$RNA@counts)

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
                   npcs = 50,
                   verbose = FALSE)
  
  seuObj <- Seurat::RunTSNE(object = seuObj,
                    dims = 1:50,
                    check_duplicates = FALSE,
                    seed.use = 1,
                    perplexity = 30,
                    verbose = FALSE)  
  return(seuObj)
  
}

# identify cells that are low quality
identify_low_quality_cells <- function(seuObj) {
  # take max deviance per cell - maximum of singlet and doublet model deviance
  singlet_dev <- seuObj@meta.data$singlet_dev
  doublet_dev <- seuObj@meta.data$singlet_dev + seuObj@meta.data$doublet_dev_imp
  max_dev <- pmax(singlet_dev, doublet_dev, na.rm=T)
  max_dev[is.na(max_dev)] <- max(max_dev, na.rm=T)
  
  # identify bimodal distribution/subset of cells with lower max deviance
  d <- density(max_dev)
  low_quality_threshold <- optimize(approxfun(d$x,d$y),interval=c(0,1))$minimum


  # plot to check that threshold is reasonable
  hist_data_set <- as.data.frame(pmax(singlet_dev, doublet_dev))
  colnames(hist_data_set) <- "max_deviance"  

  max_deviance_plot <- ggplot2::ggplot(hist_data_set,  ggplot2::aes(x=max_deviance)) +
    ggplot2::geom_histogram(data=subset(hist_data_set,max_deviance < low_quality_threshold), col='darkred', fill="firebrick2") +
    ggplot2::geom_histogram(data=subset(hist_data_set,max_deviance >= low_quality_threshold), col='dodgerblue4', fill="dodgerblue3") +
    ggplot2::geom_vline(xintercept=low_quality_threshold, linetype='dashed') +
    ggplot2::ggtitle("low quality cells")

  print(max_deviance_plot)
  
  low_quality_cells <- seuObj@meta.data$barcode[which(max_dev < low_quality_threshold)]
  
  cat(sprintf("identified %d low quality cells", length(low_quality_cells)))
  
  return(low_quality_cells)

}

# Identify doublets
classify_doublets <- function(seuObj, low_quality_cells) {
  # use singlet deviance, doublet deviance improvement, and fraction of genes detected in a given cell
  # to identify cells that are doublets
  doublet_factors <- cbind.data.frame(seuObj@meta.data$singlet_dev, seuObj@meta.data$doublet_dev_imp, seuObj@meta.data$cell_det_rate)
  rownames(doublet_factors) <- seuObj@meta.data$barcode
  colnames(doublet_factors) <- c("singlet_dev", "doublet_dev_imp", "cell_det_genes")
  # log10 transform the doublet deviance imrpvement
  doublet_factors$doublet_dev_imp <- log10(doublet_factors$doublet_dev_imp + 0.001)
  # set NA values to 0
  doublet_factors[is.na(doublet_factors)] <- 0
  
  # remove low quality cells before classifying cells as singlet or doublet
  if(length(low_quality_cells)>0) {
    lq_inds <- which(rownames(doublet_factors) %in% low_quality_cells)
    doublet_factors <- doublet_factors[-lq_inds,]
  }
  
  # use a gaussian mixture model to separate doublets and singlets
  mres <- mclust::Mclust(data = doublet_factors, G = 2)
  dq_clust <- which.max(mres$parameters$mean['doublet_dev_imp', ])
  mres_class <- data.frame(doublet_prob = mres$z[,dq_clust])
  doublets <- rownames(doublet_factors)[which(mres_class$doublet_prob>.75)]
  
  doublet_factors$doublet_probability <- mres_class$doublet_prob
  doublet_plot <-  ggplot2::ggplot(doublet_factors,  ggplot2::aes(singlet_dev, doublet_dev_imp, color=doublet_probability)) +
    ggplot2::geom_point() +
    ggplot2::ggtitle("doublet probability")
  
  plot(mres, what='classification')
  print(doublet_plot)
  
  return(doublets)
}

# add the low quality and doublet classifications for each cell
# stored as 'cell_quality' in the meta data
add_cell_quality_classifications <- function(seuObj, low_quality_cells, doublets, max_mito_perc = 0.2) {
  seuObj@meta.data$cell_quality <- "normal"
  doublet_ind <- which(seuObj@meta.data$barcode %in% doublets)
  seuObj@meta.data[doublet_ind,]$cell_quality <- "doublet"
  if(length(low_quality_cells)>0) {
    low_quality_ind <- which(seuObj@meta.data$barcode %in% low_quality_cells)
    seuObj@meta.data[low_quality_ind,]$cell_quality <- "low_quality"
  }
  
  # also classify cells with high mitochondrial percentage as low quality
  seuObj@meta.data$cell_quality[which(seuObj@meta.data$percent.mito > max_mito_perc)] <- "low_quality"
  
  cell_quality_classifications <-  ggplot2::ggplot(seuObj@meta.data,   ggplot2::aes(singlet_dev, doublet_dev_imp, color=cell_quality)) +
    ggplot2::geom_point() +  ggplot2::ggtitle("cell quality classifications")
  
  print(cell_quality_classifications)
  
  return(seuObj)
}