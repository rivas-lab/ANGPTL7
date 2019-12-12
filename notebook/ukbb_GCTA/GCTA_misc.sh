#!/bin/bash
set -beEuo pipefail

# [[ ${GRM_WB:-} -eq 1 ]]  && return || readonly GRM_WB='/oak/stanford/groups/mrivas/ukbb24983/cal/grm/ukb24983_cal_cALL_v2_hg19.white_british'
[[ ${GRM_CAL_BFILE:-}  -eq 1 ]]  && return || readonly GRM_CAL_BFILE="/scratch/groups/mrivas/ukbb/24983/cal/pgen/ukb24983_cal_cALL_v2_hg19"
[[ ${GRM_WB:-} -eq 1 ]]          && return || readonly GRM_WB='/scratch/groups/mrivas/ukbb/24983/cal/grm/ukb24983_cal_cALL_v2_hg19.white_british'
[[ ${DATA_D:-} -eq 1 ]]          && return || readonly DATA_D='/oak/stanford/groups/mrivas/projects/ANGPTL7/ukbb_GCTA'
[[ ${GCTA_phe:-} -eq 1 ]]        && return || readonly GCTA_phe="${DATA_D}/IOP_glaucoma.without.colnames.phe"
[[ ${GCTA_phe_cols:-} -eq 1 ]]   && return || readonly GCTA_phe_cols="${DATA_D}/IOP_glaucoma.phe.colnames.txt"
# [[ ${GRM_WB_ANGPTL7:-} -eq 1 ]]  && return || readonly GRM_WB="/oak/stanford/groups/mrivas/projects/ANGPTL7/ukbb_GCTA/ukb24983_cal_cALL_v2_hg19.${pop}.${var_list}"
[[ ${GRM_WB_ANGPTL7:-} -eq 1 ]]  && return || readonly GRM_WB_ANGPTL7="$PI_SCRATCH/projects/ANGPTL7/ukbb_GCTA/ukb24983_cal_cALL_v2_hg19.white_british.rs200058074_rs28991002_rs28991009_rs143435072"

get_keep () {
    local pop=$1
    echo "/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe"
}

gcta_grm_part_extract () {
    local bfile=$1
    local var_list=$2
    local out=$3
    local pop=$4
    local threads=$5
    local batch_num=$6
    local n_parts=$7

    local keep=$(get_keep $pop)

    echo ${var_list} \
    | tr "_" "\n" \
    | gcta64 \
        --extract /dev/stdin \
        --make-grm-part ${n_parts} ${batch_num} \
        --bfile ${bfile} \
        --keep ${keep} \
        --thread-num ${threads} \
        --out ${out}
}

gcta_grm_combine () {
    local part=$1
    local n_parts=$2
    local out=$3
    
    for ext in log grm.id grm.bin grm.N.bin ; do
        for i in $(seq -w $n_parts) ; do
            file_part="${part}.part_${n_parts}_${i}.${ext}"
            if [ -f ${file_part} ] ; then
                cat ${file_part}
            else
                echo [warning] missing_file: ${file_part} >&2
            fi
        done > ${out}.${ext}
    done
}

get_phe_name () {
    local pheno_idx=$1
    cat ${GCTA_phe_cols} | awk -v idx=${pheno_idx} '(NR == idx + 2){print $0}'
}

show_idx_combination (){
    local n=$1
    for i in $(seq 1 $n) ; do for j in $(seq $i $n) ; do echo $i $j ; done ; done | awk '$1 != $2'
}

gcta_HEreg () {
    local grm=$1
    local pheno_file=$2
    local pheno_idx=$3
    local out=$4
    local threads=$5

    gcta64 --grm ${grm} --pheno ${pheno_file} --mpheno ${pheno_idx} --HEreg --out ${out} --thread-num ${threads}
}

gcta_HEreg_bivar () {
    local grm=$1
    local pheno_file=$2
    local pheno_idx1=$3
    local pheno_idx2=$4
    local out=$5
    local threads=$6

    gcta64 --grm ${grm} --pheno ${pheno_file} --HEreg-bivar ${pheno_idx1} ${pheno_idx2} --out ${out} --thread-num ${threads}
}

get_HEreg_file () {
    local file_prefix=$1
    echo "/oak/stanford/groups/mrivas/projects/ANGPTL7/ukbb_GCTA/white_british/${file_prefix}.HEreg"
}

gcta_HEreg_show () {
    local file=$1
    # this function dumps "V(G)/Vp" estimates, SE, and its p-value.
    # 1. HE-CP V(G)/Vp Estimate
    # 2. HE-CP V(G)/Vp SE_Jackknife
    # 3. HE-CP V(G)/Vp P_Jackknife
    # 4. HE-SD V(G)/Vp Estimate
    # 5. HE-SD V(G)/Vp SE_Jackknife
    # 6. HE-SD V(G)/Vp P_Jackknife
    cat $file | grep 'V(G)/Vp' | awk -v OFS='\t' '{print $2, $4, $5}' | tr "\n" "\t"
}

get_HEreg_bivar_file () {
    local file_prefix1=$1
    local file_prefix2=$2
    echo "/oak/stanford/groups/mrivas/projects/ANGPTL7/ukbb_GCTA/white_british/${file_prefix1}_${file_prefix2}.HEreg"
}

gcta_HEreg_bivar_show () {
    local file=$1
    # this function dumps "V(G)/Vp" estimates, SE, and its p-value.
    # 1. HE-CP rG Estimate
    # 2. HE-CP rG SE_Jackknife
    # 3. HE-CP rG SE_OLS
    
    cat $file | grep 'rG' | awk -v OFS='\t' '{print $2, $4, $3}'
}

# phe
#. /oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/INI2005254.phe