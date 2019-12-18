suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# tab_file <- '/oak/stanford/groups/mrivas/private_data/ukbb/16698/release1/download/ukb9444.csv'
# out_file <- 'glaucoma_relevant_raw_data_ukb9444.tsv'

tab_file <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/9797/download/ukb9797.tab'
out_file <- 'glaucoma_relevant_raw_data_ukb9797.tsv'

icd_keep <- c(
    'H40',
    'H41',
    'H400',
    'H401',
    'H402',
    'H403',
    'H404',
    'H405',
    'H406',
    'H408',
    'H409',
    'H42',
    'H420',    
    'H428',
    'Q150'
)

# select_cols <- c(
#     'eid',
#     paste0('20002-0.', 0:28),
#     paste0('20002-1.', 0:28),
#     paste0('20002-2.', 0:28),
#     paste0('41202-0.', 0:379)
# )

select_cols <- c(
    'f.eid',
    paste0('f.20002.0.', 0:28),
    paste0('f.20002.1.', 0:28),
    paste0('f.20002.2.', 0:28),
    paste0('f.41202.0.', 0:379)
)

fread(tab_file, select=select_cols) %>%
rename('IID' = 'f.eid') %>% 
# rename('IID' = 'eid') %>% 
gather(field, val, -IID) %>% 
drop_na(val) %>%
# mutate(field = str_replace_all(field, '-', '.')) %>%
mutate(field = str_replace_all(field, '^f.', '')) %>%
separate(
    field, c('UKB_Field_ID', 'UKB_time_idx', 'UKB_array_idx')
) %>%
filter(
    ((UKB_Field_ID == '20002') & (val == '1277')) |
    ((UKB_Field_ID == '41202') & (val %in% icd_keep))
) %>% 
select(-UKB_time_idx, -UKB_array_idx) %>% 
unique() %>%
mutate(UKB_Field_ID_val = paste0('F', UKB_Field_ID, '_', val)) %>% 
select(-UKB_Field_ID, -val) %>%
mutate(value = TRUE) %>%
spread(UKB_Field_ID_val, value, fill=FALSE) %>%
fwrite(out_file, sep='\t')
