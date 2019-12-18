suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

read_col_names <- function(tab_file){
    df <- fread(cmd=paste0(
        "head -n1 ", tab_file,
        " | tr '\t' '\n' ",
        " | sed -e 's/^f.//g' ",
        " | tr '.' '\t' ",
        " | awk 'NF==3'"
    ))
    colnames(df) <- c('UKB_Field_ID', 'UKB_time_idx', 'UKB_array_idx')
    df
}

get_select_cols <- function(tab_file, UKB_Field_IDs){
    tab_file %>% read_col_names() %>% 
    filter(
        UKB_Field_ID %in% UKB_Field_IDs
    ) %>%
    mutate(
        colname = paste('f', UKB_Field_ID, UKB_time_idx, UKB_array_idx, sep='.')
    ) %>%
    select(colname) %>% 
    pull()
}

raw_to_wide <- function(raw_df, UKB_Field_IDs, keep_vals){
    raw_df %>%
    rename('IID' = 'f.eid') %>% 
    gather(field, val, -IID) %>% 
    drop_na(val) %>%
    mutate(field = str_replace_all(field, '^f.', '')) %>%
    separate(
        field, c('UKB_Field_ID', 'UKB_time_idx', 'UKB_array_idx')
    ) %>%
    filter(
        ((UKB_Field_ID %in% UKB_Field_IDs[['self_reported']]) & (val == keep_vals[['self_reported']])) |
        ((UKB_Field_ID %in% UKB_Field_IDs[['ICD']]) & (val %in% keep_vals[['ICD']]))
    ) %>% 
    select(-UKB_time_idx, -UKB_array_idx) %>% 
    unique() %>%
    mutate(UKB_Field_ID_val = paste0('F', UKB_Field_ID, '_', val)) %>% 
    select(-UKB_Field_ID, -val) %>%
    mutate(value = TRUE) %>%
    spread(UKB_Field_ID_val, value, fill=FALSE)    
}

tab_file <- '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/9796/9797/download/ukb9797.tab'
out_file <- 'ukb9797.glaucoma.tsv'

UKB_Field_IDs <- list(
    self_reported = c('20002'),
    ICD  = c('41202','42104','40001','41078','41143','41142','41104','41079','41105','40002','40006','41204')
)

keep_vals <- list(
    self_reported = c('1277'),
    ICD = c('H40', 'H41', 'H400', 'H401', 'H402', 'H403', 'H404', 'H405', 'H406', 'H408', 'H409', 'H42', 'H420', 'H428', 'Q150')
)

select_cols <- c('f.eid', get_select_cols(tab_file, unname(unlist(UKB_Field_IDs))))
print(length(select_cols))

wide_df <- tab_file %>% 
fread(select=select_cols) %>%
raw_to_wide(UKB_Field_IDs, keep_vals)

wide_df %>% fwrite(out_file, sep='\t')
