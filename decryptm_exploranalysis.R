library(tidyverse)

experimental_summary <- read_tsv("experiment_summary.tsv")

collectri_tfs <- read_tsv('collectri.tsv') %>% pull(source) %>% unique()

# fitted data
fitted_params <- read_csv('fit_params_peptide_rep_MD.csv')
# 1192 psites, 892 with fit, 142 r2>0.8 (54 prots)


results_tfs_longf <- read_csv('detected_peptides_long_MD.csv') %>% dplyr::rename('rep' = rep, 'Drug' = drug, 'Cell_line' = cell_line) %>% filter(Gene %in% collectri_tfs)
concentrations <- results_tfs_longf %>% distinct(conc) %>% arrange(conc) %>% pull(conc)
# full_data <- full_join(fitted_curves, results_tfs_longf, by = c('Gene', 'Drug', 'Cell_line'))

logistic_function <- function(x, ec50, slope, top, bottom){
    return((top - bottom) / (1 + 10 ** (slope * (x - ec50))) + bottom)
}

# create curves from the logistic function and the fitted parameters
fitted_curves <- fitted_params %>% 
    filter(rsq>0.8) %>%
    group_by(Index, Gene, Drug, Cell_line, rep) %>%
    group_split() %>%
    purrr::map(., function(x){
        data <- tibble(x_fit = seq(0, 10000, 1))
        data$y_fit <- logistic_function(log10(data$x_fit), unique(x$log_ec50), unique(x$slope), unique(x$top), unique(x$bottom))
        data$Index <- unique(x$Index)
        data$Gene <- unique(x$Gene)
        data$Drug <- unique(x$Drug)
        data$Cell_line <- unique(x$Cell_line)
        data$rep <- unique(x$rep)
        return(data)
    }) %>% bind_rows() %>% filter(!is.na(y_fit))

# open a pdf devvice, then plot the curves (four curves corresponding to the four replicates per plot). Also plot the real measurements as points. Add r2 tag to each plot

fitted_params <- fitted_params %>% ungroup()
fitted_curves <- fitted_curves %>% ungroup()




unique_combinations <- fitted_params %>% distinct(Index, Gene, Drug, rep, rsq) %>% arrange(desc(rsq))

combs_over_threshold <- unique_combinations %>% filter(rsq >= 0.8)
combs_under_threshold <- unique_combinations %>% filter(rsq < 0.8) %>% slice_sample(n=40)

# iterate over all combinations of gene, drug, peptide and rep
pdf('fitted_curves.pdf', width = 10, height = 7)
combs_over_threshold %>%
  group_by(Index, Gene, Drug, rep) %>%
  group_walk(~{
    # Aquí puedes acceder a cada grupo con .x y .y
    # .x es el dataframe del grupo actual
    # .y es un dataframe con las claves del grupo actual
    gene <- .y$Gene
    drug <- .y$Drug
    index <- .y$Index
    rep <- .y$rep

    p <- fitted_curves %>% 
        filter(Gene == gene, Drug == drug, Index == index, rep==!!rep) %>% 
        ggplot(aes(x = x_fit, y = y_fit)) +
        geom_line() +
        scale_x_log10() +
        ggtitle(paste0(index, ' ', drug, ' ', rep)) +
        cowplot::theme_cowplot() +
        # plot r2
        annotate("text", x = 0.1, y = 1, label = paste0('r2 = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Index == index, rep == !!rep) %>% pull(rsq)), 2))) +
        annotate("text", x = 0.1, y = 1.1, label = paste0('top = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Index == index, rep == !!rep) %>% pull(top)), 2))) +
        annotate("text", x = 0.1, y = 1.2, label = paste0('log ec50 = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Index == index, rep == !!rep) %>% pull(log_ec50)), 2))) +
        annotate("text", x = 0.1, y = 1.3, label = paste0('bottom = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Index == index, rep == !!rep) %>% pull(bottom)), 2))) +
        annotate("text", x = 0.1, y = 1.4, label = paste0('slope = ', round(unique(fitted_params %>% filter(Gene == gene, Drug == drug, Index == index, rep == !!rep) %>% pull(slope)), 2))) +
        theme(legend.position = 'none')
    
    # print real measurements in the same plot
    p2 <- p + geom_point(data = results_tfs_longf %>% filter(Gene == gene, Drug == drug, Index == index, rep==!!rep), aes(x = conc, y = ratio), alpha = 0.5)

    print(p2)
  })


combs_under_threshold %>%
    group_by(Index, Gene, Drug, rep) %>%
    group_walk(~{
        gene <- .y$Gene
        drug <- .y$Drug
        index <- .y$Index
        rep <- .y$rep
        
        p <- results_tfs_longf %>% 
            filter(Gene == gene, Drug == drug, Index == index, rep==!!rep) %>% 
            ggplot(aes(x = conc, y = ratio)) +
            # add title
            ggtitle(paste0(index, ' ', drug, ' ', rep)) +
            geom_point() +
            scale_x_log10() +
            cowplot::theme_cowplot() +
            theme(legend.position = 'none')
        print(p)
    })
dev.off()


# 33 genes w peptides with at least r2 > 0.8 and log_ec50 < 2
formatted_targets <- fitted_params %>% filter(log_ec50<2, rsq>0.8) %>% 
    mutate(bio_context = paste(Cell_line, Drug, sep = '_')) %>%
    arrange(Gene, Drug, Cell_line) %>%
    mutate(direction = sign(bottom-top)) %>%
    distinct(Gene, Drug, Cell_line, direction, .keep_all = TRUE) %>%
    arrange(log_ec50) %>%
    relocate(direction, .after = log_ec50) %>%
    group_by(Gene, Drug, Cell_line) %>%
    dplyr::filter(log_ec50 == min(log_ec50), .bygroup=TRUE) %>%
    ungroup() %>%
    select(Gene, bio_context, direction) %>%
    mutate(direction = as.numeric(direction)) %>%
    pivot_wider(id_cols = bio_context, names_from=Gene, values_from = direction) %>%
    column_to_rownames('bio_context')

non_phosphorylated_tfs <- collectri_tfs[!collectri_tfs %in% names(formatted_targets)] %>%
tibble(TF = ., A431_Afatinib=NA, A431_Gefitinib = NA, A431_Dasatinib = NA) %>%
column_to_rownames('TF') %>% t() %>% data.frame()

formatted_targets <- merge(formatted_targets, non_phosphorylated_tfs, by='row.names') %>% column_to_rownames('Row.names')

write.table(formatted_targets, 'proteomics_targets.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)




