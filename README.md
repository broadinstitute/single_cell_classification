# single_cell_classification
Methods to use SNPs or gene expression to classify single cell RNAseq to reference profiles

These methods were used in the paper: 'Multiplexed single-cell profiling of post-perturbation transcriptional responses to define cancer vulnerabilities and therapeutic mechanism of action'

## Running single_cell_classification:

This method was run on macOS High Sierra v10.13.6 using RStudio version 1.2.5033 and R version 3.6.2 (2019-12-12).

* clone or download the git repo
* install the R package dependencies using the install_packages.R script
* open single_cell_classification as your project root directory
* To run classifications use source(here::here('src', 'run_SNP_classification.R')); run_SNP_classification()
* To run the QC methods use source(here::here('src', 'run_QC.R')); run_all_QC()

run_SNP_classification() and run_all_QC() read in data from the folder 'data' within the single_cell_classification folder. There is a test data set of 588 cells originating from 5 cell lines. The test data includes the cell barcodes, gene list, expression matrix, and SNP ref and alt allele count matrices for the single cell data and the SNP ref and alt allele count matrices for the 5 reference profiles. On the test data provided run_SNP_classification() and run_all_QC() each take <1 min to run. 

## Files needed:

These scripts are intended to run on 10X single cell RNAseq processed using Cell Ranger (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). 

For run_SNP_classification, these files are needed:

* bulk_ref.csv
* bulk_alt.csv
* sc_ref.csv
* sc_alt.csv
* barcodes.tsv

Matrices of reference and alternate allele counts at given SNP sites are needed for reference cell lines and single cell data. The methods assume that the alt and ref matrices are samples x SNP sites. For bulk RNAseq, freebayes (https://github.com/ekg/freebayes) was run with forced calling at a set of 100,000 SNP sites, then we combined the freebayes vcf file for each sample using load_data_helpers.R/combine_bulk_reference_profiles to create ref and alt allele count matrices. For the single cell RNAseq, ref and alt allele matrices were produced using the method scAlleleCount (https://github.com/barkasn/scAlleleCount). The single cell barcodes, as output from the cellranger count method, are also required.


To run the QC scripts or gene_expression_classification these files, as output from the cellranger count method, are needed:

* matrix.mtx
* genes.tsv
* barcodes.tsv
* classifications.csv - which is the output of run_SNP_classification()

## Output: 

run_SNP_classification, using defaults, will output a matrix containing the reference sample classifications for each cell, as well as other classification metrics, such as:

* singlet_ID : the most likely reference sample classification for that cell NA,
* singlet_dev : fraction of deviance explained by the top reference sample
* singlet_dev_z : (zscored) fraction of deviance explained by the top reference sample
* singlet_margin : difference between the fraction of deviance explained by the most likely reference sample and the second most likely reference sample
* singlet_z_margin :  difference between the (zscored) fraction of deviance explained by the most likely reference sample and the second most likely reference sample
* doublet_z_margin :  difference between the (zscored) fraction of deviance explained by the second most likely reference sample and the third most likely reference sample
* doublet_dev_imp : the difference between the fraction of deviance explained by the doublet model and the fraction of deviance explained by the singlet model (measure of whether it is a doublet or singlet)
* doublet_CL1 : the most likely reference cell line if this cell is a doublet
* doublet_CL2 : the next most likely reference cell line if this cell is a doublet
* tot_reads : total reads in this cell at the SNP sites used for classification 
* num_SNPs : number of SNPs detected in this cell of the SNP sites used for classification

run_QC outputs a Seurat object and includes cell quality classifications for each cell (in the meta.data object), classifying each cell as:

* normal : cells that are used for downstream analyis
* doublet : cell is more likely a multiplet, discared for downstream analysis
* low quality : cell is low quality in terms of RNAseq quality metrics or in terms of our ability to classify it as one of the reference samples, discarded for downstream analysis
* empty droplet : cells, with distinct gene expression profiles, and SNP profiles that did not match to any reference cell line (or pairwise combination of cell lines) in particular, but rather resembled more a mixture of SNPs from all the in-pool cell lines, suggesting these are empty droplets containing ambient mRNA in the pool, discarded for downstream analysis
* low confidence : cells that could not be confidently classified as any of the reference samples, discarded for downstream analysis

gene_expression_classification can be used for comparison to SNP based classifications, but is not recommended for primary classification.

