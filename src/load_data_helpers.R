
# read in bulk reference and alternate allele counts used to classify single cell data
# reference profiles were created by running freebayes
load_reference_files <- function(bulk_ref_file, bulk_alt_file, pseudo_cnt = 1) {
  bulk_alt <- read.csv(bulk_alt_file, check.names = F) %>% as.data.frame()
  bulk_ref <- read.csv(bulk_ref_file, check.names = F) %>% as.data.frame()
  
  # Add dummy SNP to prevent errors with glmnet:
  bulk_alt <- rbind(bulk_alt, rep(1,ncol(bulk_alt))) %>% as.matrix()
  bulk_ref <- rbind(bulk_ref, rep(1,ncol(bulk_ref))) %>% as.matrix()
  
  # add pseudocount
  bulk_alt_p <- (bulk_alt + pseudo_cnt)/(bulk_alt + bulk_ref + 2*pseudo_cnt) %>% as.matrix()
  
  return(bulk_alt_p)
  
}

# read in single cell reference and alternate allele counts
# files are output from snpclust
load_single_cell_files <- function(sc_ref_file, sc_alt_file) {
  sc_ref <- Matrix::readMM(sc_ref_file) %>% as.matrix() %>% t()
  sc_alt <- Matrix::readMM(sc_alt_file) %>% as.matrix() %>% t()
  
  # Add dummy SNP to prevent errors with glmnet:
  sc_ref <- rbind(sc_ref, rep(1,ncol(sc_ref))) %>% as.matrix()
  sc_alt <- rbind(sc_alt, rep(1,ncol(sc_alt))) %>% as.matrix()
  
  return(list(sc_ref=sc_ref, sc_alt=sc_alt))
}

# load counts matrix, genes, and cell identities for a given experiment
load_sc_exp <- function(data_folder) {
  rMat <- Matrix::readMM(file.path(data_folder, 'matrix.mtx')) %>% as.matrix()
  genes <- readr::read_tsv(file.path(data_folder, 'genes.tsv'), col_names = F) %>% 
    set_colnames(c('Ensembl_ID', 'Gene_Symbol'))
  barcodes <- readr::read_tsv(file.path(data_folder, 'barcodes.tsv'), col_names = F) 
  colnames(rMat) <- barcodes$X1
  rownames(rMat) <- genes$Ensembl_ID
  num_CLs <- ncol(rMat)
  num_genes <- nrow(rMat)
  classification <- barcodes
  if (file.exists(file.path(data_folder, 'classifications.csv'))) {
    classification <- readr::read_csv(file.path(data_folder, 'classifications.csv')) %>% as.data.frame()
  }
  cat(sprintf('Loaded matrix with %d CLs and %d genes\n', num_CLs, num_genes))
  
  return(list(counts = rMat, gene_info = genes, cell_info = classification))
}

# convert gene identifiers from ensembl IDs to hugo symbols
convert_to_hugo_symbols <- function (dat) {
  ensemble_to_hugo <- dat$gene_info$Gene_Symbol %>% magrittr::set_names(dat$gene_info$Ensembl_ID)
  rownames(dat$counts) %<>% ensemble_to_hugo[.] %>% magrittr::set_names(NULL)
  dat$counts <- dat$counts[!is.na(rownames(dat$counts)), ]
  dat$counts <- dat$counts[!duplicated(rownames(dat$counts)), ]
  return(dat)
}

# downsample single cell ref and alt allele counts to test SNP classification
downsample_single_cell <- function(sc_alt, sc_ref, downsample_rate) {
  sc_ref <- DropletUtils::downsampleMatrix(sc_ref, downsample_rate)
  sc_alt <- DropletUtils::downsampleMatrix(sc_alt, downsample_rate)
  
  return(list(sc_ref=sc_ref, sc_alt=sc_alt))
}

# combine bulk reference profiles output from freebayes
# to make reference profile matrices
combine_bulk_reference_profiles <- function(reference_folder, snp_file) {
  # get list of reference files
  reference_files <- list.files(path = reference_folder, full.names = TRUE)
  reference_names <- gsub(".vcf", "", list.files(reference_folder))
  
  # read in SNP set used for classification
  snps <- data.table::fread(snp_file) %>% as.data.frame()
  snp_names <- paste0(snps[,1], "_", snps[,2], "_", snps[,3], "_", snps[,4])
  
  ref_count <- matrix(data=0, nrow=nrow(snp_file), ncol=length(reference_files))
  alt_count <- matrix(data=0, nrow=nrow(snp_file), ncol=length(reference_files))
  
  colnames(ref_count) <- reference_names
  rownames(ref_count) <- snp_names
  colnames(alt_count) <- reference_names
  rownames(alt_count) <- snp_names
  
  # combine freebayes output into matrices
  for(i in 1:length(reference_files)) {
    cur_bulk <- vcfR::read.vcfR(bulk_files[i])
    bulk_snps <- paste0(cur_bulk@fix[,1], "_", cur_bulk@fix[,2], "_", cur_bulk@fix[,4], "_", cur_bulk@fix[,5])
    ro_ind <- grep("RO", strsplit(cur_bulk@gt[1,1], ":")[[1]])
    ao_ind <- grep("AO", strsplit(cur_bulk@gt[1,1], ":")[[1]])
    #SNP sites with multiple base pairs
    rep_snps <- grep(",", bulk_snps)
    # get snps at these sites
    multi_snps <- data.table::tstrsplit(bulk_snps[rep_snps], ",")
    # number of alt alleles
    num_alt <- length(multi_snps)
    # get first snp
    snp1 <-  multi_snps[[1]]
    # second base
    snp_base2 <-  multi_snps[[2]]
    # third base
    if(num_alt == 3) {
      snp_base3 <- multi_snps[[3]]
    }
    # base of snp name w/out alt allele
    base_snp <- substr(snp1,1,nchar(snp1)-1)
    # snp name with second alt allele
    snp2 <- paste0(base_snp, snp_base2)
    # snps w/ 3 alt alleles
    if(num_alt == 3) {
      base3_snps <- which(!is.na(snp_base3))
      # snp names with third alt allele
      snp3 <- paste0(base_snp[base3_snps], snp_base3[base3_snps])
    }
    
    # replace comma name with first alt allele name
    bulk_snps[rep_snps] <- snp1
    # concatenate other alt alleles to end of SNP list
    bulk_snps <- c(bulk_snps, snp2)
    if(num_alt == 3) {
      bulk_snps <- c(bulk_snps, snp3)
      
    }
    
    # get alt counts
    AO <- data.table::tstrsplit(cur_bulk@gt[,2], ":", keep=ao_ind)[[1]]
    
    #SNP sites NA RO and AO
    NA_AO <- which(is.na(cur_bulk@gt[,2]) == T)
    
    # SNP sites with NA RO/AO and not multi alt alleles
    NA_AO_1 <- setdiff(NA_AO, rep_snps)
    
    # SNP sites with NA RO/AO and  multi alt alleles
    NA_rep_alleles <- which(rep_snps %in% NA_AO)
    NA_AO_3 <- c()
    if(num_alt == 3) {
      NA_three_ind <- NA_rep_alleles[which(is.na(multi_snps[[3]][NA_rep_alleles])==F)]
      NA_AO_3<- rep_snps[NA_three_ind]
    }
    NA_AO_2 <- setdiff(intersect(rep_snps, NA_AO), NA_AO_3)
    
    cur_bulk@gt[NA_AO_1,2] <- '0/0:0:0:0:0:0:0:0'
    cur_bulk@gt[NA_AO_2,2] <- '0/0:0:0:0:0:0,0:0:0'
    cur_bulk@gt[NA_AO_3,2] <- '0/0:0:0:0:0:0,0,0:0:0'
    
    #SNP sites with multi alt counts
    AO <- data.table::tstrsplit(cur_bulk@gt[,2], ":", keep=ao_ind)[[1]]
    rep_AO <- grep(",", AO)
    multi_AO <- data.table::tstrsplit(AO[rep_AO], ",")
    
    # Split multi alt counts
    AO1 <-  multi_AO[[1]]
    AO2 <-  multi_AO[[2]]
    if(num_alt == 3) {
      AO3 <- multi_AO[[3]]
      AO3 <- AO3[which(!is.na(AO3))]
    }
    AO[rep_AO] <- AO1
    AO <- c(AO, AO2)
    if(num_alt == 3) {
      AO <- c(AO, AO3)
    }
    
    AO <- as.numeric(AO)    
    AO[which(is.na(AO))] <- 0
    
    # get ref count
    RO <- as.numeric(data.table::tstrsplit(cur_bulk@gt[,2], ":", keep=ro_ind)[[1]])
    RO[which(is.na(RO))] <- 0
    
    # concatenate same RO for snps with second and third alt alleles
    RO <- c(RO, RO[rep_snps])
    if(num_alt == 3) {
      RO <- c(RO, RO[rep_snps[base3_snps]])
    }
    
    # fill in count matrix
    used_inds <- which(bulk_snps %in% rownames(ref_count))
    ref_count[bulk_snps[used_inds],i] <- RO[used_inds]
    alt_count[bulk_snps[used_inds],i] <- AO[used_inds]
  }
}
