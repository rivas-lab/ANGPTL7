suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

df <- "ukb9797.glaucoma.tsv" %>% fread()
phe_df <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/highconfidenceqc/phe/HC276.phe' %>% fread()
names(phe_df) <- c('FID', 'IID', 'phe')

setequal(
    intersect(df$IID, phe_df$IID), 
    phe_df %>% filter(phe == 2) %>% select(IID) %>% pull() 
) %>% print()
