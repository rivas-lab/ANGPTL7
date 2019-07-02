
```
tabix -h /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/9796/24611/white_british/ukb24983_v2_hg19.INI2734.genotyped.PHENO1.glm.linear.gz 1:11252357-11253688 | awk -v OFS='\t' -v sep=':' '{print $1 sep $2 sep $4 sep $5, $2, $8, $9, $12}' | sed -e "s%#CHROM:POS:REF:ALT%#ID%g" > INI2734.ANGPTL7.tsv
tabix -h /oak/stanford/groups/mrivas/ukbb24983/cal/gwas/10136/21731/white_british/ukb24983_v2_hg19.INI2405.genotyped.PHENO1.glm.linear.gz 1:11252357-11253688 | awk -v OFS='\t' -v sep=':' '{print $1 sep $2 sep $4 sep $5, $2, $8, $9, $12}' | sed -e "s%#CHROM:POS:REF:ALT%#ID%g" > INI2405.ANGPTL7.tsv
```

