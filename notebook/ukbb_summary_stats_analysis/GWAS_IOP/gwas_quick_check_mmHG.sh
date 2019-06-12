#!/bin/bash
set -beEuo pipefail

#cat /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/10136/21731/phe/INI5255.phe | awk '($3 == -9){print $1}' | sort -b | join -1 1 -2 1 /dev/stdin <(cat /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe | cut -f1 | sort -b) | awk '{print $0, $0}' > ukb24983_white_british_no_IOP.phe


#var_id=$1

#var_id='rs28991009'

#echo $var_id \
cat ../../ukbb_individual_level_data_analysis/ANGPTL7.protein-altering.vars.lst \
| plink2 \
  --bpfile /oak/stanford/groups/mrivas/private_data/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2_hg19 \
  --chr 1-22 \
  --covar /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_GWAS_covar.phe \
  --covar-name age sex Array PC1-PC4 \
  --glm genotypic firth-fallback hide-covar omit-ref \
  --memory 25600 \
  --pheno /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/phe/HC276.phe \
  --remove /oak/stanford/groups/mrivas/ukbb24983/sqc/ukb24983_v2.not_used_in_pca.phe \
  --threads 4 \
  --extract /dev/stdin \
  --keep /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe \
  --out ukb24983_v2_hg19.INI5255.genotyped.mmHg

#  --pheno-quantile-normalize \
#  --keep ukb24983_white_british_no_IOP.phe \
#  --extract /oak/stanford/groups/mrivas/ukbb24983/sqc/both_array_variants.txt \
#  --out /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/extras/highconfidenceqc/white_british/ukb24983_v2_hg19.HC276.genotyped.both_arrays

