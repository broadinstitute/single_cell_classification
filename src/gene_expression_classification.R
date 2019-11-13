# Gene expression classifications
# intended for comparing to SNP based classifications, but not intended as primary classifications

source(here::here('src', 'load_data_helpers.R'))

#PARAMS
n_highvar_genes_compare <- 5000 #number of genes to use for comparing CCLE and cluster avgs
vtr <- c()
param_range <- 10

# can be downloaded from here: https://depmap.org/portal/download/
# rows are samples and columns are genes of the format 'Hugo symbol (Ensemble ID)'
CCLE_GE <- read_csv(here::here('data', 'CCLE_expression_full.csv')) %>% column_to_rownames('X1') %>% as.matrix()

dat <- load_sc_exp(here::here('data'))
seuObj <- create_seurat_obj(dat, use_symbols = FALSE)


n_cls <- nlevels(seuObj)
n_pcs <- 2*n_cls
if(n_cls > 25-param_range & n_cls < 25+param_range) {
  clust_res <- 1
} else if(n_cls > 50-param_range & n_cls < 50+param_range) {
  clust_res <- 2
} else if(n_cls > 100-param_range & n_cls < 100+param_range) {
  clust_res <- 4
} else {
  stop('Not implemented')
}

#use just good singlets
cq <- Seurat::FetchData(seuObj, vars = c('cell_quality'))
seuObj <- seuObj[, which(cq$cell_quality == 'normal')]
seuObj <- Seurat::ScaleData(object = seuObj, vars.to.regress = vtr)
seuObj <- Seurat::FindVariableFeatures(object = seuObj,
                                       nfeatures = 1e5,
                                       selection.method = 'vst')
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

#cluster cells
seuObj <- Seurat::FindNeighbors(seuObj, reduction = 'pca',
                                dims = 1:n_pcs,
                                k.param = 10, 
                                force.recalc = TRUE,
                                verbose = FALSE)
seuObj <- Seurat::FindClusters(seuObj, resolution = clust_res, verbose = FALSE)

#compute summed counts across cells for each cluster, then CPM transform
clust_cpm <- plyr::laply(levels(seuObj), function(clust) {
  Matrix::rowSums(Seurat::GetAssayData(seuObj, slot = 'counts')[, names(Seurat::Idents(seuObj)[Seurat::Idents(seuObj) == clust])])
}) %>%
  t() %>%
  edgeR::cpm(prior.count = 1, log = TRUE)

in_pool_lines <- unique(seuObj@meta.data$singlet_ID)
in_pool_lines_depmap <- unique(seuObj@meta.data$DepMap_ID)

colnames(CCLE_GE) <- stringr::str_match(colnames(CCLE_GE), '\\((ENSG[0-9]+)\\)')[,2]
sc_ensemble <- seuObj@misc$Ensembl_ID
common_genes <- intersect(colnames(CCLE_GE), sc_ensemble)

high_var_genes <- apply(CCLE_GE[in_pool_lines_depmap, match(common_genes, colnames(CCLE_GE))], 2, sd) %>% 
  sort(decreasing = T) %>% 
  head(5000) %>% 
  names()  

#mean-center both CCLE and cluster expression profiles
CCLE_GE_norm <- CCLE_GE[in_pool_lines_depmap, high_var_genes] %>% scale(center = T, scale = F) %>% t() 
clust_cpm_norm <- clust_cpm[match(high_var_genes, seuObj@misc$Ensembl_ID),] %>% t() %>% scale(center = T, scale = F) %>% t()

cc <- cor(CCLE_GE_norm, clust_cpm_norm)
best_matches <- in_pool_lines_depmap[apply(cc, 2, which.max)]
GE_class <- best_matches[as.numeric(Seurat::Idents(seuObj))]

perc_correct <- mean(GE_class == seuObj@meta.data$DepMap_ID)

seuObj@meta.data$GE_class <- GE_class

seuObj@meta.data %<>% mutate(agree = GE_class == DepMap_ID)

df <- cbind.data.frame(seuObj@reductions$tsne@cell.embeddings, seuObj@meta.data)

ggplot2::ggplot(df, ggplot2::aes(tSNE_1, tSNE_2, color = agree)) + 
  ggplot2::geom_point() +
  ggplot2::ggtitle("Agreement between SNP and gene expression classifications")
