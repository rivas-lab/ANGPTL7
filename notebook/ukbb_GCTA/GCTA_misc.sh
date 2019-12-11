#!/bin/bash
set -beEuo pipefail

# [[ ${GRM_WB:-} -eq 1 ]]  && return || readonly GRM_WB='/oak/stanford/groups/mrivas/ukbb24983/cal/grm/ukb24983_cal_cALL_v2_hg19.white_british'
[[ ${GRM_WB:-} -eq 1 ]]         && return || readonly GRM_WB='/scratch/groups/mrivas/ukbb/24983/cal/grm/ukb24983_cal_cALL_v2_hg19.white_british'
[[ ${DATA_D:-} -eq 1 ]]         && return || readonly DATA_D='/oak/stanford/groups/mrivas/projects/ANGPTL7/ukbb_GCTA'
[[ ${GCTA_phe:-} -eq 1 ]]       && return || readonly GCTA_phe="${DATA_D}/IOP_glaucoma.without.colnames.phe"
[[ ${GCTA_phe_cols:-} -eq 1 ]]  && return || readonly GCTA_phe_cols="${DATA_D}/IOP_glaucoma.phe.colnames.txt"

get_phe_name () {
    local pheno_idx=$1
    cat ${GCTA_phe_cols} | awk -v idx=${pheno_idx} '(NR == idx + 2){print $0}'
}

show_idx_combination (){
    local n=$1
    for i in $(seq 1 $n) ; do for j in $(seq $i $n) ; do echo $i $j ; done ; done | awk '$1 != $2'
}

gcta_HEreg () {
    local pheno_file=$1
    local pheno_idx=$2
    local out=$3
    local threads=$4

    gcta64 --grm ${GRM_WB} --pheno ${pheno_file} --mpheno ${pheno_idx} --HEreg --out ${out} --thread-num ${threads}
}

gcta_HEreg_bivar () {
    local pheno_file=$1
    local pheno_idx1=$2
    local pheno_idx2=$3
    local out=$4
    local threads=$5

    gcta64 --grm ${GRM_WB} --pheno ${pheno_file} --HEreg-bivar ${pheno_idx1} ${pheno_idx2} --out ${out} --thread-num ${threads}
}


# phe
#. /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/INI2005254.phe