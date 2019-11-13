# single_cell_classification
Methods to use SNPs or gene expression to classify single cell RNAseq to reference profiles

These methods were used in the paper: 'Pan-cancer single cell RNA-seq uncovers recurring programs of cellular heterogeneity'

## Files needed:

These scripts are intended to run on 10X single cell RNAseq processed using Cell Ranger (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger). 

For run_SNP_classification, these files are needed:

* bulk_ref.csv
* bulk_alt.csv
* sc_ref.csv
* sc_alt.csv
* barcodes.tsv

Matrices of reference and alternate allele counts at given SNP sites are needed for reference cell lines and single cell data. The methods assume alt and ref matrices for the reference samples are SNP sites x reference samples and alt and ref matrices for the single cell data are cells x SNP sites. For bulk RNAseq, freebayes (https://github.com/ekg/freebayes) was run with forced calling at a set of 100,000 SNP sites, then we combined the freebayes vcf file for each sample using load_data_helpers.R/combine_bulk_reference_profiles to create ref and alt allele count matrices. For the single cell RNAseq, ref and alt allele matrices were produced using the method snpclust (https://github.com/10XGenomics/single-cell-3prime-snp-clustering). The single cell barcodes, as output from the cellranger count method, are also required.


To run the QC scripts or gene_expression_classification these files, as output from the cellranger count method, are needed:

* matrix.mtx
* genes.tsv
* barcodes.tsv

## Output: 

run_SNP_classification, using defaults, will output a matrix containing the reference sample classifications for each cell, as well as other classification metrics, such as:

* singlet_ID : the most likely reference sample classification for that cell
* singlet_dev : fraction of deviance explained by the top reference sample
* singlet_dev_z : (zscored) fraction of deviance explained by the top reference sample
* singlet_z_margin :  difference between the (zscored) fraction of deviance explained by the most likely reference sample and the second most likely reference sample
* doublet_dev_imp : the difference between the fraction of deviance explained by the doublet model and the fraction of deviance explained by the singlet model (measure of whether it is a doublet or singlet)
* doublet_CL1 : the most likely reference cell line if this cell is a doublet
* doublet_CL2 : the next most likely reference cell line if this cell is a doublet
* tot_reads : total reads in this cell at the SNP sites used for classification 
* num_SNPs : number of SNPs detected in this cell of the SNP sites used for classification

run_QC outputs cell quality classifications for each cell, classifying each cell as:

* normal : cells that are used for downstream analyis
* doublet : cell is more likely a multiplet, discarded for downstream analysis
* low quality : cell is low quality in terms of RNAseq quality metrics or in terms of our ability to classify it as one of the reference samples, discarded for downstream analysis

gene_expression_classification can be used for comparison to SNP based classifications, but is not recommended for primary classification. It produces similar output to SNP based classification. 
