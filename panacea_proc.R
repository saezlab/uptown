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
  



  
  targets_drug <- panacea_targets %>% filter(cmpd == treatment) %>% select(target) %>% pull()

  panacea_treatment_filtered$regulation <- ifelse(panacea_treatment_filtered$logFC > 1.5 & panacea_treatment_filtered$padj < 0.05, "Upregulated",
                                    ifelse(panacea_treatment_filtered$logFC < -1.5 & panacea_treatment_filtered$padj < 0.05, "Downregulated", "No Change"))

  panacea_treatment_filtered$target <- panacea_treatment_filtered$gene_symbol %in% targets_drug & (panacea_treatment_filtered$regulation == "Upregulated" | panacea_treatment_filtered$regulation == "Downregulated")

  volcano_plot <- ggplot(panacea_treatment_filtered, aes(x = logFC, y = -log10(padj), color = regulation)) +
    geom_point(alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "blue") +
    geom_vline(xintercept = -1.5, linetype = "dashed", color = "blue") +
    geom_text(data = subset(panacea_treatment_filtered, target & (regulation == "Upregulated" | regulation == "Downregulated")),
              aes(label = gene_symbol), vjust = -0.5, hjust = 0.5) +
    theme_minimal() +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No Change" = "black")) +
    xlab("Log2 Fold Change") +
    ylab("-Log10 P-value") +
    theme(legend.position = "none") +
    facet_grid(cell_line~treatment)
  
  plot_list[[i]] <- volcano_plot
}

library(egg)
volcanos <- egg::ggarrange(plots = plot_list, ncol = 32)
ggsave('panacea_volcanos_tmm_limma.png', volcanos, width = 100, height = 40, dpi = 200, limitsize=FALSE)


full_targets <- panacea_targets %>% select(cmpd, target, rank)

des_targets <- full_targets %>% filter(rank == 1)

sources_tsv <- des_targets %>% select(-rank) %>% mutate(sign = -1)

write_tsv(sources_tsv, 'panacea_sources.tsv')

offtargets_df <- full_targets %>% filter(rank != 1) %>% select(-rank)
offtargets_df

write_tsv(offtargets_df, 'panacea_offtargets.tsv')


panacea_targets_df <- tibble()
for(treatment in treatments){
  panacea_treatment_filtered <- panacea_de_filtered %>%
    select(gene_symbol, contains(treatment)) %>%
    pivot_longer(cols = -gene_symbol, names_to = 'pipeline', values_to = 'value') %>%
    separate_wider_delim(pipeline, delim = '__', names = c('statparam', 'status', 'pipeline', 'biocontext', 'main_dataset', 'subset')) %>%
    separate_wider_delim(biocontext, delim = '_', names = c('cell_line', 'treatment', 'v', 'dmso1', 'dmso2')) %>%
    select(-dmso1, -v, -dmso2, -status, -pipeline, -main_dataset, -subset) %>%
    pivot_wider(names_from = statparam, values_from = value) %>%
    filter(logFC > 1.5 & padj < 0.05) %>%
    arrange(desc(abs(stat)))
  
  panacea_targets_df <- bind_rows(panacea_targets_df, panacea_treatment_filtered)

}

write_tsv(panacea_targets_df, 'panacea_targets.tsv')



