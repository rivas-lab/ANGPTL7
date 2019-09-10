#!/bin/bash
set -beEuo pipefail

for GBE_ID in INI5254 INI5262 INI5263 ; do
    bash ~/repos/rivas-lab/ukbb-tools/07_LDSC/scripts/ukb_ldsc_munge_rg.sh \
    data/ukb24983_v2_hg19.INI5255.genotyped.PHENO1.glm.linear.gz \
    data/ukb24983_v2_hg19.${GBE_ID}.genotyped.PHENO1.glm.linear.gz \
    ANGPTL7
done
