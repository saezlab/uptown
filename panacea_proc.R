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
  select(ID, contains('__filtered__tmm+limma')) %>%
  dplyr::rename(gene_symbol = ID) # %>%
  # dplyr::rename(logFC = contains('logFC'),
  #               pval = contains('padj'),
  #               stat = contains('stat'))
plot_list <- list() 

for (i in 1:length(treatments)){
  treatment <- treatments[i]
  panacea_treatment_filtered <- panacea_de_filtered %>%
    select(gene_symbol, contains(treatment)) %>%
    pivot_longer(cols = -gene_symbol, names_to = 'pipeline', values_to = 'value') %>%
    separate_wider_delim(pipeline, delim = '__', names = c('statparam', 'status', 'pipeline', 'biocontext', 'main_dataset', 'subset')) %>%
    separate_wider_delim(biocontext, delim = '_', names = c('cell_line', 'treatment', 'v', 'dmso1', 'dmso2')) %>%
    select(-dmso1, -v, -dmso2, -status, -pipeline, -main_dataset, -subset) %>%
    pivot_wider(names_from = statparam, values_from = value)

  
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

panacea_treatment_filtered <- panacea_de_filtered %>%
  select(gene_symbol, contains('AFATINIB')) %>%
  pivot_longer(cols = -gene_symbol, names_to = 'pipeline', values_to = 'value') %>%
  separate_wider_delim(pipeline, delim = '__', names = c('statparam', 'status', 'pipeline', 'biocontext', 'main_dataset', 'subset')) %>%
  separate_wider_delim(biocontext, delim = '_', names = c('cell_line', 'treatment', 'v', 'dmso1', 'dmso2')) %>%
  select(-dmso1, -v, -dmso2) %>%
  pivot_wider(names_from = statparam, values_from = value)

targets_drug <- panacea_targets %>% filter(cmpd == 'AFATINIB') %>% select(target) %>% pull()

panacea_treatment_filtered$regulation <- ifelse(panacea_treatment_filtered$logFC > 1.5 & panacea_treatment_filtered$padj < 0.05 & panacea_treatment_filtered$gene_symbol %in% targets_drug, "Upregulated",
                                    ifelse(panacea_treatment_filtered$logFC < -1.5 & panacea_treatment_filtered$padj < 0.05 & panacea_treatment_filtered$gene_symbol %in% targets_drug, "Downregulated", "No Change"))

ggplot(panacea_treatment_filtered, aes(x = logFC, y = -log10(padj), color = regulation)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "blue") +
  theme_minimal() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No Change" = "black")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  facet_wrap(~cell_line)



panacea_targets_in_genes <- c()
for (target in panacea_targets_list){
  if (target %in% panacea_genes){
    panacea_targets_in_genes <- c(panacea_targets_in_genes, target)
  }
}

# write a for loop wrapping an if else statement assessing if a string is in a list:

for (target in panacea_targets_list){
  if (target %in% panacea_genes){
    panacea_targets_in_genes <- c(panacea_targets_in_genes, target)
  }
}



# Create a new variable to indicate upregulation or downregulation
panacea_de_filtered$regulation <- ifelse(panacea_de_filtered$logFC > 1.5 & panacea_de_filtered$pval < 0.05, "Upregulated",
                                    ifelse(panacea_de_filtered$logFC < -1.5 & panacea_de_filtered$pval < 0.05, "Downregulated", "No Change"))

# Create the volcano plot
volcano_plot <- ggplot(panacea_de_filtered, aes(x = logFC, y = -log10(pval), color = regulation)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "blue") +
  geom_vline(xintercept = -1.5, linetype = "dashed", color = "blue") +
  theme_minimal() +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No Change" = "black")) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value")

# View the plot
print(volcano_plot)
