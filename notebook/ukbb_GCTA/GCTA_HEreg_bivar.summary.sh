#!/bin/bash
set -beEuo pipefail

# source "$(dirname ${SRCNAME})/GCTA_misc.sh"
source "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ANGPTL7/notebook/ukbb_GCTA/GCTA_misc.sh"
phe_num=7 # the total number of phenotypes

echo "#GBE_ID1 GBE_ID2 rG SE_Jackknife SE_OLS" | tr " " "\t"

for i in $(seq $(perl -e "print($phe_num * ($phe_num - 1) / 2)")) ; do
    pheno_idx1="$(show_idx_combination $phe_num | awk -v nr=${i} '(NR == nr){print $1}')"
    pheno_idx2="$(show_idx_combination $phe_num | awk -v nr=${i} '(NR == nr){print $2}')"
    GBE_ID1="$(get_phe_name ${pheno_idx1})"
    GBE_ID2="$(get_phe_name ${pheno_idx2})"
    HEreg_file=$(get_HEreg_bivar_file ${GBE_ID1} ${GBE_ID2})
    gcta_HEreg_bivar_show ${HEreg_file} | awk -v gbe1=$GBE_ID1 -v gbe2=$GBE_ID2 -v OFS='\t' '{print gbe1, gbe2, $0}'
done
