#!/bin/bash
set -beEuo pipefail

for pop in white_british non_british_white african e_asian s_asian ; do
    comm -2 -3  \
        <(cat /oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_${pop}.phe | tr "\t" "_"  | sort ) \
        <(cat phe/INI*.phe | grep -v 'FID' | awk -v OFS='_' '{print $1, $2}' | sort -u) \
        | tr "_" "\t" > phe/ukb24983_${pop}_noIOP.keep
done

wc phe/ukb24983_*_noIOP.keep

