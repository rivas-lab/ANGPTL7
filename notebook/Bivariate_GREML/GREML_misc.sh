#!/bin/bash
set -beEuo pipefail

[[ ${GRM_WB:-} -eq 1 ]]  && return || readonly GRM_WB='/oak/stanford/groups/mrivas/ukbb24983/cal/grm/ukb24983_cal_cALL_v2_hg19.white_british'

gcta_reml () {
    pheno=$1
    out=$2
    threads=$3
    gcta64 --grm $GRM_WB --pheno $pheno --reml --out $out --thread-num ${threads}
}



# phe
#. /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/INI2005254.phe