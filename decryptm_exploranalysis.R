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
fitted_params <- read_csv('fit_params_peptide_rep_nc.csv')


results_tfs_longf <- read_csv('detected_peptides_long.csv') %>% dplyr::rename('rep' = rep, 'Drug' = drug, 'Cell_line' = cell_line, 'Peptide' = `Modified Peptide`) %>% filter(Gene %in% collectri_tfs)
concentrations <- results_tfs_longf %>% distinct(conc) %>% arrange(conc) %>% pull(conc)
# full_data <- full_join(fitted_curves, results_tfs_longf, by = c('Gene', 'Drug', 'Cell_line'))

logistic_function <- function(x, ec50, slope, top, bottom){
    return((top - bottom) / (1 + 10 ** (slope * (x - ec50))) + bottom)
}

# create curves from the logistic function and the fitted parameters
fitted_curves <- fitted_params %>% 
    filter(rsq>0.8) %>%
    group_by(Gene, Drug, Cell_line, rep, Peptide) %>%
    group_split() %>%
    purrr::map(., function(x){
        data <- tibble(x_fit = seq(0, 10000, 1))
        data$y_fit <- logistic_function(data$x_fit, unique(x$log_ec50), unique(x$slope), unique(x$top), unique(x$bottom))
        data$Gene <- unique(x$Gene)
        data$Drug <- unique(x$Drug)
        data$Cell_line <- unique(x$Cell_line)
        data$rep <- unique(x$rep)
        data$Peptide <- unique(x$Peptide)
        return(data)
    }) %>% bind_rows() %>% filter(!is.na(y_fit))

# open a pdf devvice, then plot the curves (four curves corresponding to the four replicates per plot). Also plot the real measurements as points. Add r2 tag to each plot

fitted_params <- fitted_params %>% ungroup()
fitted_curves <- fitted_curves %>% ungroup()




unique_combinations <- fitted_curves %>% distinct(Gene, Drug, Peptide, rep) %>% arrange(Gene, Drug, Peptide, rep)

# iterate over all combinations of gene, drug, peptide and rep
pdf('fitted_curves.pdf', width = 10, height = 7)
unique_combinations %>%
  group_by(Gene, Drug, Peptide, rep) %>%
  group_walk(~{
    # Aquí puedes acceder a cada grupo con .x y .y
    # .x es el dataframe del grupo actual
    # .y es un dataframe con las claves del grupo actual
    gene <- .y$Gene
    drug <- .y$Drug
    peptide <- .y$Peptide
    rep <- .y$rep
    if ((fitted_curves %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep==!!rep) %>% nrow() == 0) & (results_tfs_longf %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep==!!rep) %>% nrow() > 0)){
        p <- results_tfs_longf %>% 
            filter(Gene == gene, Drug == drug, Peptide == peptide, rep==rep) %>% 
            ggplot(aes(x = conc, y = ratio)) +
            # add title
            ggtitle(paste0(gene, ' ', drug, ' ', peptide, ' ', rep)) +
            geom_point() +
            scale_x_log10() +
            cowplot::theme_cowplot() +
            theme(legend.position = 'none')
        print(p)
    } else if (results_tfs_longf %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep==!!rep) %>% nrow() == 0){
        next
    } else{
    p <- fitted_curves %>% 
        filter(Gene == gene, Drug == drug, Peptide == peptide, rep==!!rep) %>% 
        ggplot(aes(x = x_fit, y = y_fit)) +
        geom_line() +
        scale_x_log10() +
        ggtitle(paste0(gene, ' ', drug, ' ', peptide, ' ', rep)) +
        cowplot::theme_cowplot() +
        # plot r2
        annotate("text", x = 0.1, y = 1, label = paste0('r2 = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep == !!rep) %>% pull(rsq)), 2))) +
        annotate("text", x = 0.1, y = 1.1, label = paste0('top = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep == !!rep) %>% pull(top)), 2))) +
        annotate("text", x = 0.1, y = 1.2, label = paste0('log ec50 = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep == !!rep) %>% pull(log_ec50)), 2))) +
        annotate("text", x = 0.1, y = 1.3, label = paste0('bottom = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep == !!rep) %>% pull(bottom)), 2))) +
        annotate("text", x = 0.1, y = 1.4, label = paste0('slope = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep == !!rep) %>% pull(slope)), 2))) +
        theme(legend.position = 'none')
    
    # print real measurements in the same plot
    p2 <- p + geom_point(data = results_tfs_longf %>% filter(Gene == gene, Drug == drug, Peptide == peptide, rep==!!rep), aes(x = conc, y = ratio), alpha = 0.5)
    print(p2)
    }
  })
dev.off()


# create a subset of 100 random combinations of gene, drug, peptide and rep which were not fitted
unique_combinations <- fitted_curves %>% distinct(Gene, Drug, Peptide, rep) %>% arrange(Gene, Drug, Peptide, rep)



