Prerequisites:

python 2.7.13 or 3.6.1
numpy at least 1.16.4
pandas at least 0.24.2
scipy at least 1.2.2
colorama
rpy2 compatible with the version of python (2.8.6 for python 2.7 and 2.9.2 for python 3.6)

Usage: gene_based_test.py [-h] --file FILE_PATH --pop POP --pheno PHENO
                          --metadata_path METADATA_PATH
                          [--R_var {independent,similar} [{independent,similar} ...]]
                          [--variants {pcv,pav,ptv} [{pcv,pav,ptv} ...]]
                          [--maf_thresh MAF_THRESHES [MAF_THRESHES ...]]
                          [--out_folder OUT_FOLDER]

Argument Details:
  -h, --help            show this help message and exit
  --file FILE_PATH      path to tab-separated file containing summary statistics.
                               
                                 format of necessary columns:
                                
                                 #CHROM   POS   REF  ALT  BETA/OR  SE
                                 
                                 
  --pop POP             name of the population. used in output file naming.
                                 
  --pheno PHENO         name of the phenotype. used in output file naming.
                                 
  --metadata_path METADATA_PATH
                        path to tab-separated file containing:
                                 variants,
                                 gene symbols,
                                 consequences,
                                 MAFs,
                                 and LD independence info.
                                 
                                 such a file is provided in this repository.
                               
                                 format:
                                 
                                 V       gene_symbol     most_severe_consequence maf  ld_indep
                                 1:69081:G:C     OR4F5   5_prime_UTR_variant     0.000189471     False
                                
  --R_var {independent,similar} [{independent,similar} ...]
                        type(s) of model across variants. 
                                 options: independent, similar (default: independent). can run both.
                                 independent model is akin to the dispersion test; similar model is akin to burden.

  --variants {pcv,pav,ptv} [{pcv,pav,ptv} ...]
                        variant set(s) to consider. 
                                 options: proximal coding [pcv], 
                                          protein-altering [pav], 
                                          protein truncating [ptv] 
                                          (default: ptv). can run multiple.

  --maf_thresh MAF_THRESHES [MAF_THRESHES ...]
                        which MAF threshold(s) to use. must be valid floats between 0 and 1 
                                 (default: 0.01).

  --out_folder OUT_FOLDER
                        folder to which output(s) will be written (default: current folder).
                                 if folder does not exist, it will be created.

Example usage (used to generate files provided):
python gene_based_test.py --file /oak/stanford/groups/mrivas/projects/ANGPTL7/ukbb_gwas/white_british/ukb24983_v2_hg19.INI5255.genotyped.glm.linear.gz --pop white_british --pheno INI5255 --metadata_path ukb_cal-consequence_wb_maf_gene_ld_indep.tsv --R_var independent similar --variants pav ptv
