library(tidyverse)

experimental_summary <- read_tsv("experiment_summary.tsv")

collectri_tfs <- read_tsv('collectri.tsv') %>% pull(source) %>% unique()

# read network and get list of nodes
network <- read_tsv('network_collectri.sif')
nodes <- network %>% 
    select('# source', target)

# append sources and targets in a list
nodes_list <- list(nodes$`# source`, nodes$target) %>% unlist() %>% unique()

fitted_params <- read_csv('fit_params_peptide_all_bound.csv')
# fitted data
ec50_info <- fitted_params %>% 
    filter(rsq>=0.8) %>%
    group_by(Gene, Drug, Cell_line) %>%
    summarise(log_ec50 = min(log_ec50)) %>%
    filter(Gene %in% nodes_list) %>%
    filter(!Gene %in% collectri_tfs) #summarised at protein level, best phosphorylation site per protein
# summarise log_ec50 values per gene, drug and cell line

write_tsv(ec50_info, 'ec50_info.tsv')

# number of collectri tfs with rsq > 0.8
fitted_params %>% 
    # filter(rsq>=0.8) %>%
    filter(Gene %in% nodes_list) %>%
    filter(Gene %in% collectri_tfs) %>%
    group_by(Gene, Drug, Cell_line)

results_tfs_longf <- read_csv('detected_peptides_long_MD.csv') %>% dplyr::rename('rep' = rep, 'Drug' = drug, 'Cell_line' = cell_line) %>% filter(Gene %in% collectri_tfs)
concentrations <- results_tfs_longf %>% distinct(conc) %>% arrange(conc) %>% pull(conc)
# full_data <- full_join(fitted_curves, results_tfs_longf, by = c('Gene', 'Drug', 'Cell_line'))

logistic_function <- function(x, ec50, slope, top, bottom){
    return((top - bottom) / (1 + 10 ** (slope * (x - ec50))) + bottom)
}

# create curves from the logistic function and the fitted parameters
fitted_curves <- fitted_params %>% 
    filter(rsq>0.8) %>%
    filter(Gene %in% collectri_tfs) %>%
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




unique_combinations <- fitted_params %>% distinct(Index, Gene, Drug, rep, rsq) %>% arrange(desc(rsq)) %>% filter(Gene %in% collectri_tfs)

combs_over_threshold <- unique_combinations %>% filter(rsq >= 0.8)
combs_under_threshold <- unique_combinations %>% filter(rsq < 0.8) %>% slice_sample(n=100)

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


phosphorylated_prots <- fitted_params %>% pull(Gene) %>% unique()

formatted_targets <- fitted_params %>% 
    filter(rsq>0.8) %>% 
    filter(Gene %in% collectri_tfs) %>%
    filter(Gene %in% nodes_list) %>%
    mutate(bio_context = paste(Cell_line, Drug, sep = '_')) %>%
    group_by(Index, Gene, Drug, Cell_line, bio_context) %>%
    summarise(log_ec50 = mean(log_ec50), top = min(top), bottom=min(bottom)) %>%
    mutate(direction = sign(bottom-top)) %>%
    distinct(Index, Gene, Drug, Cell_line, direction, .keep_all = TRUE) %>%
    relocate(direction, .after = log_ec50) %>%
    group_by(Gene, Drug, Cell_line, bio_context) %>%
    summarise(log_ec50 = min(log_ec50), direction = min(direction)) %>%
    filter(direction != 0) %>%
    ungroup() %>%
    group_by(Drug, Cell_line, bio_context) %>%
    arrange(log_ec50, .by_group=TRUE) %>%
    group_by(bio_context) %>%
    select(Gene, bio_context, direction) %>%
    mutate(direction = as.numeric(direction)) %>%
    pivot_wider(id_cols = bio_context, names_from=Gene, values_from = direction) %>%
    column_to_rownames('bio_context')

# plot distribution of log_ec50, facet per biocontext
fitted_params %>% 
    filter(rsq>0.8) %>% 
    filter(Gene %in% collectri_tfs) %>%
    mutate(bio_context = paste(Cell_line, Drug, sep = '_')) %>%
    group_by(Index, Gene, Drug, Cell_line, bio_context) %>%
    summarise(log_ec50 = mean(log_ec50), top = min(top), bottom=min(bottom)) %>%
    mutate(direction = sign(bottom-top)) %>%
    distinct(Index, Gene, Drug, Cell_line, direction, .keep_all = TRUE) %>%
    relocate(direction, .after = log_ec50) %>%
    group_by(Gene, Drug, Cell_line, bio_context) %>%
    summarise(log_ec50 = min(log_ec50), direction = min(direction)) %>%
    filter(direction != 0) %>%
    ungroup() %>%
    group_by(Drug, Cell_line, bio_context) %>%
    arrange(log_ec50, .by_group=TRUE) %>%
    group_by(bio_context) %>%
    slice_min(n=25, order_by=log_ec50)   %>%
    ggplot(aes(x=log_ec50)) +
    geom_histogram(bins=30) +
    facet_wrap(~bio_context, ncol=3) +
    cowplot::theme_cowplot() +
    theme(legend.position = 'none')

non_changing_tfs <- phosphorylated_prots[phosphorylated_prots %in% collectri_tfs & !phosphorylated_prots %in% names(formatted_targets)] %>%
tibble(TF = ., A431_Afatinib=NA, A431_Gefitinib = NA, A431_Dasatinib = NA) %>%
column_to_rownames('TF') %>% t() %>% data.frame()

formatted_targets <- merge(formatted_targets, non_changing_tfs, by='row.names') %>% column_to_rownames('Row.names')

write.table(formatted_targets, 'proteomics_targets.tsv', sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)



