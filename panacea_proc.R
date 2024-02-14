library(tidyverse)

panacea_de <- read_tsv('../../flop/flop_results/diffexp/GSE186341__deresults.tsv')
panacea_targets <- read_csv('./PANACEA/panacea_gold_standard.csv') %>%
  group_by(cmpd) %>%
  arrange(rank, .by_group = TRUE)


treatments <- panacea_targets %>% select(cmpd) %>% distinct() %>% pull()

cell_lines <- panacea_de %>% select(contains('BAFETINIB')) %>% names(.) %>%
  as_tibble() %>% separate_wider_delim(value, delim = '__', names_sep = 'a') %>%
  select(valuea4) %>%
  separate_wider_delim(valuea4, delim = '_', names_sep = 'a') %>%
  distinct(valuea4a1) %>% pull()

panacea_de_filtered <- panacea_de %>%
  select(ID, contains('__filtered__NA+deseq2')) %>%
  dplyr::rename(gene_symbol = ID) # %>%
  # dplyr::rename(logFC = contains('logFC'),
  #               pval = contains('padj'),
  #               stat = contains('stat'))
plot_list <- list()
mega_df <- tibble()

for (i in 1:length(treatments)){
  treatment <- treatments[i]
  panacea_treatment_filtered <- panacea_de_filtered %>%
    select(gene_symbol, contains(treatment)) %>%
    pivot_longer(cols = -gene_symbol, names_to = 'pipeline', values_to = 'value') %>%
    separate_wider_delim(pipeline, delim = '__', names = c('statparam', 'status', 'pipeline', 'biocontext', 'main_dataset', 'subset')) %>%
    separate_wider_delim(biocontext, delim = '_', names = c('cell_line', 'treatment', 'v', 'dmso1', 'dmso2')) %>%
    select(-dmso1, -v, -dmso2, -status, -pipeline, -main_dataset, -subset) %>%
    pivot_wider(names_from = statparam, values_from = value)
  
  # if there are less than 30 genes with abs(logFC) > 1.5 and pvalue<0.05, per cell line and treatment, remove the cell line-treatment combo
  sig_cell_lines <- panacea_treatment_filtered %>%
    group_by(cell_line, treatment) %>%
    filter(abs(logFC) > 1.5 & padj < 0.05) %>%
    count() %>% filter (n>30) %>%
    pull(cell_line) %>% unique()

  sig_data <- panacea_treatment_filtered %>% filter(cell_line %in% sig_cell_lines)
  
  mega_df <- bind_rows(mega_df, sig_data)
}

write_tsv(mega_df, 'panacea_de.tsv')
  

full_targets <- panacea_targets %>% select(cmpd, target, rank)

des_targets <- full_targets %>% filter(rank == 1)

sources_tsv <- des_targets %>% select(-rank) %>% mutate(sign = -1)

write_tsv(sources_tsv, 'panacea_sources.tsv')

offtargets_df <- full_targets %>% filter(rank != 1) %>% select(-rank)
offtargets_df

write_tsv(offtargets_df, 'panacea_offtargets.tsv')



