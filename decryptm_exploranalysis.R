library(tidyverse)

experimental_summary <- read_tsv("experiment_summary.tsv")

results_full <- read_tsv('decryptm/all_phospho_processed_full_results/tmt-report/abundance_protein_MD.tsv')
results_tfs <- read_tsv('decryptm/all_phospho_processed_tfs_results/tmt-report/abundance_protein_MD.tsv')
collectri_tfs <- read_tsv('collectri.tsv') %>% pull(source) %>% unique()

sample_ids <- tibble(sample_id = c(paste0('sample-0', seq(1,9,1)), 'sample-10', 'sample-11'))

concentrations <- experimental_summary %>% select(starts_with('Conc')) %>% distinct() %>% t() %>% as.data.frame() %>% rownames_to_column('conc_id') %>% as_tibble() %>% rename(conc = V1) %>% bind_cols(sample_ids)

results_full_longf <- results_full %>% 
    pivot_longer(cols = starts_with('sample-'), names_to = 'sample_id', values_to = 'ratio') %>% 
    left_join(concentrations, by = 'sample_id')  %>%
    filter(Gene %in% collectri_tfs)

write_tsv(results_full_longf, 'phosphorliated_prots_clean.tsv')

results_tfs_longf <- results_tfs %>% 
    pivot_longer(cols = starts_with('sample-'), names_to = 'sample_id', values_to = 'abundance') %>% 
    left_join(concentrations, by = 'sample_id') %>%
    filter(Gene %in% collectri_tfs)

detected_proteins_full <- results_full_longf %>% distinct(Gene) %>% pull()

detected_proteins_tfs <- results_tfs_longf %>% distinct(Gene) %>% pull()

# plot share of detected proteins present in collectri
#71 proteins detected when using the full proteome
detected_proteins_full %>% 
    intersect(collectri_tfs) %>% 
    length()

#71 proteins detected when using the tfs
detected_proteins_tfs %>% 
    intersect(collectri_tfs) %>% 
    length()

# plot the activity (y axis) vs concentration (x axis) for each protein
# plot the activity (y axis) vs concentration (x axis) for each protein
results_tfs_longf %>% 
    filter(Gene %in% collectri_tfs) %>% 
    ggplot(aes(x = conc, y = abundance, color = Gene)) +
    geom_line() +
    scale_x_log10() +
    cowplot::theme_cowplot() +
    theme(legend.position = 'none') +
    facet_wrap(~Gene, ncol = 8)



# fitted data
fitted_data <- read_tsv('phosphorylated_prots_fitted.tsv')

logistic_function <- function(x, A, B, C, D){
    return(A + (C - A) / (1 + exp(-(x - B) / D)))
}

# per TF, plot the fitted curve and the raw data
fitted_data %>% 
    filter(Gene %in% collectri_tfs) %>% 
    mutate(fitted = logistic_function(conc, A, B, C, D)) %>% 
    ggplot(aes(x = conc, y = ratio, color = Gene)) +
    geom_point() +
    geom_line(aes(y = fitted),color='red') +
    scale_x_log10() +
    cowplot::theme_cowplot() +
    theme(legend.position = 'none') +
    facet_wrap(~Gene, ncol = 8)
