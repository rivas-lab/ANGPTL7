suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
library(UpSetR)

in_file <- 'ukb9797.glaucoma.tsv'

keep_files <- list(
    WB_noIOP = '/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/ukb24983_white_british_noIOP.keep',
    WB = '/oak/stanford/groups/mrivas/ukbb24983/sqc/population_stratification/ukb24983_white_british.phe' 
)
out_plot_prefix <- 'ukb9797.glaucoma'

read_keep <- function(file){
    df <- file %>% fread(sep='\t')
    colnames(df) <- c('FID', 'IID')
    df
}

wide_to_plot <- function(wide_df){
    wide_df %>%
    na_if(FALSE) %>%
    gather(field, value, -IID) %>% 
    drop_na(value) %>%
    mutate(field = str_replace_all(field, '^F', '')) %>%
    separate(
        field, c('UKB_Field_ID', 'coded_value')
    ) %>%
    mutate(coded_value = str_replace_all(coded_value, '1277', 'Glaucoma')) %>%
    select(-UKB_Field_ID, -value) %>%
    # select(-UKB_Field_ID) %>%    
    unique()
    #  %>%
    # spread(coded_value, value, fill=FALSE)    
}

get_IIDs_for_upset_plot <- function(plot_df, val){
    plot_df %>% filter(coded_value == val) %>% select(IID) %>% pull()
}

upset_plot_wrapper <- function(plot_df){
    upset(
        fromList(list(
            'Self-reported' = plot_df %>% get_IIDs_for_upset_plot('Glaucoma'),
            # 'ICD-10 H40'    = plot_df %>% get_IIDs_for_upset_plot('H40'),
            # 'ICD-10 H41'    = plot_df %>% get_IIDs_for_upset_plot('H41'),
            'ICD-10 H400'   = plot_df %>% get_IIDs_for_upset_plot('H400'),
            'ICD-10 H401'   = plot_df %>% get_IIDs_for_upset_plot('H401'),
            'ICD-10 H402'   = plot_df %>% get_IIDs_for_upset_plot('H402'),
            'ICD-10 H403'   = plot_df %>% get_IIDs_for_upset_plot('H403'),
            'ICD-10 H404'   = plot_df %>% get_IIDs_for_upset_plot('H404'),
            'ICD-10 H405'   = plot_df %>% get_IIDs_for_upset_plot('H405'),
            'ICD-10 H406'   = plot_df %>% get_IIDs_for_upset_plot('H406'),
            'ICD-10 H408'   = plot_df %>% get_IIDs_for_upset_plot('H408'),
            'ICD-10 H409'   = plot_df %>% get_IIDs_for_upset_plot('H409'),
            # 'ICD-10 H42'    = plot_df %>% get_IIDs_for_upset_plot('H42'),
            'ICD-10 H420'   = plot_df %>% get_IIDs_for_upset_plot('H420'),
            'ICD-10 H428'   = plot_df %>% get_IIDs_for_upset_plot('H428'),
            'ICD-10 Q150'   = plot_df %>% get_IIDs_for_upset_plot('HQ50')
        )), 
        order.by = "freq",
        nsets = 15, nintersects = 500,
        number.angles = 0, 
        point.size = 2, line.size = .5, 
        mainbar.y.label = "Number of case individuals", 
        sets.x.label = "# cases per data source", 
        text.scale = c(1, 1, 1, 1, 1, .8)
    )    
}

df <- fread(in_file)

pdf(file=paste0(out_plot_prefix, '.all.pdf'), onefile=FALSE, height = 6, width=18)
df %>%
wide_to_plot() %>% 
upset_plot_wrapper()
dev.off()

pdf(file=paste0(out_plot_prefix, '.white_british.pdf'), onefile=FALSE, height = 6, width=18)
df %>%
filter(IID %in% read_keep(keep_files[['WB']])$IID) %>%
wide_to_plot() %>% 
upset_plot_wrapper()
dev.off()

pdf(file=paste0(out_plot_prefix, '.white_british_noIOP.pdf'), onefile=FALSE, height = 6, width=18)
df %>%
filter(IID %in% read_keep(keep_files[['WB_noIOP']])$IID) %>%
wide_to_plot() %>% 
upset_plot_wrapper()
dev.off()
