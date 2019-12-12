#!/bin/bash
set -beEuo pipefail

# source "$(dirname ${SRCNAME})/GCTA_misc.sh"
source "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ANGPTL7/notebook/ukbb_GCTA/GCTA_misc.sh"
phe_num=7 # the total number of phenotypes

echo "#GBE_ID V(G)/Vp_CP SE_Jackknife_CP P_Jackknife_CP V(G)/Vp_SD SE_Jackknife_SD P_Jackknife_SD" | tr " " "\t"

for i in $(seq $phe_num) ; do
    GBE_ID="$(get_phe_name $i)"
    HEreg_file=$(get_HEreg_file "${GBE_ID}.rs200058074_rs28991002_rs28991009_rs143435072")
    gcta_HEreg_show ${HEreg_file} | awk -v gbe=$GBE_ID -v OFS='\t' '{print gbe, $0}'
done
