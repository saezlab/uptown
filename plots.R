library(ggplot2)
library(ggridges)
library(tidyverse)
library(cowplot)
library(egg)
library(ggpattern)
panacea_data <- read_csv('panacea_graphdata_results.csv')

# change the random column so that it contains True when Random==real and False otherwise
panacea_data$Random <- ifelse(panacea_data$Random == 'real', 'False', 'True')

methods <- panacea_data %>% pull(Methods) %>% unique() %>% sort()
biocontexts <- panacea_data %>% 
  group_by(Biocontext) %>%
# filter out the biocontext for which all rows from the column offtarget_counts is na
  filter(all(!is.na(`offtarget_count`))) %>%
  filter() %>% pull(Biocontext) %>% unique() %>% sort()

panacea_data <- panacea_data %>% 
    filter(Biocontext %in% biocontexts) %>%
    mutate(Methods = factor(Methods, levels = methods))
    # keep biocontexts present in the biocontext list

biocontexts_order <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    group_by(Biocontext) %>%
    summarise(mean_biocontext_diff = mean(mean_diff, na.rm=TRUE)) %>%
    arrange(mean_biocontext_diff)%>% pull(Biocontext)

panacea_data <- panacea_data %>% 
    mutate(Biocontext = factor(Biocontext, levels = biocontexts_order))

panacea_dotplot <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    filter(!is.na(sd_diff)) %>%
    # create dotplot, x is method, y is bio_context, color is mean_diff, size is sd_diff (bigger sd smaller size)
    ggplot(aes(x = Methods, y = Biocontext, color = mean_diff, size = sd_diff)) +
    geom_point() +
    #scale color threshold: <0 is blue, +0 is red
    scale_color_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0, limits = c(-50, 50)) +
    scale_size_continuous(range = c(0, 5), trans='reverse') +
    labs(color = 'Mean % diff\nreal - random', size = 'Null sd') +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=10), 
    axis.text.y = element_text(size=10),
    axis.title = element_text(size=10),
    legend.position = 'right',
    text = element_text(size=10),
    legend.box = "horizontal")

legend <- get_legend(panacea_dotplot) %>% ggpubr::as_ggplot()

panacea_dotplot <- panacea_dotplot + theme(legend.position = 'none')

top_annotation <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    group_by(Methods) %>%
    mutate(mean_method_diff = mean(mean_diff, na.rm=TRUE)) %>%
    ggplot(aes(x = Methods, y = mean_diff, fill = mean_method_diff)) +
    # color boxplot by mean y value
    geom_boxplot() +
    scale_fill_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0, limits = c(-50, 50)) +
    geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.75) +
    theme_cowplot() +
    theme(text=element_text(size=10),
          axis.text.y = element_text(size=10),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=10),
          legend.position = 'none')

right_annotation <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    group_by(Biocontext) %>%
    mutate(mean_biocontext_diff = mean(mean_diff, na.rm=TRUE)) %>%
    ggplot(aes(y = Biocontext, x = mean_diff, fill = mean_biocontext_diff)) +
    # color boxplot by mean y value
    geom_boxplot() +
    scale_fill_gradient2(low = 'blue', mid = 'grey', high = 'red', midpoint = 0, limits = c(-50, 50)) +
    geom_vline(xintercept = 0, linetype = 'dashed', linewidth = 0.75) +
    theme_cowplot() +
    theme(text=element_text(size=10),
      axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size=10),
          axis.text.x = element_text(size=10),
          legend.position = 'none'
          )



# The main plot is in the bottom left, the top annotation is above it, and the right annotation is to the right
# The placeholder is in the top right to balance the layout
arranged_plots <- ggarrange(
  top_annotation, legend,
  panacea_dotplot, right_annotation,
  ncol = 2, nrow = 2,
  widths = c(2.5, 1), # Adjust as necessary
  heights = c(1, 4)) # Adjust as necessary


ggsave('plots/panacea_dotplot.png', plot = arranged_plots, device = 'png', width = 19, height = 29, units = 'cm', dpi = 200)

# numberss

panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    group_by(Methods) %>%
    summarise(mean_method_diff = mean(mean_diff, na.rm=TRUE)) %>%
    arrange(desc(mean_method_diff)) %>% print(n=40)

panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    group_by(Biocontext) %>%
    summarise(mean_method_diff = mean(mean_diff, na.rm=TRUE)) %>% arrange(desc(mean_method_diff)) %>% print(n=40)



offtargets_df <- read_tsv('panacea_offtargets.tsv') %>% group_by(cmpd) %>%
    summarise(count = n())
    
disagreement_df <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    summarise(mean = mean(`perc_offtarget`, na.rm=TRUE), sd = sd(`perc_offtarget`, na.rm=TRUE)) %>%
    group_by(Biocontext, Methods) %>%
    summarise(mean_diff = -diff(mean), sd_diff = sd) %>%
    filter(!is.na(sd_diff)) %>%
    separate(Biocontext, c('cell_line', 'treatment'), sep='_', remove = FALSE) %>%
    left_join(offtargets_df, by = c('treatment' = 'cmpd'))

disagreement_plot <- ggplot(disagreement_df, aes(x = count, y = mean_diff, size = sd_diff)) +  # Removed size = sd_diff from here to apply only to geom_point
  geom_point() +  # Apply size = sd_diff only to geom_point
  # geom_ribbon(data = min_max_df, aes(x = count, ymin = min_mean_diff, ymax = max_mean_diff), inherit.aes = FALSE, fill = "grey", alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 0.75) +
  scale_color_gradient2(low = 'red', mid = 'black', high = 'red', midpoint = 0, limits = c(-50, 50)) +
  theme_cowplot()

plot_bigdispersion <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    filter(Biocontext == 'H1793_OSIMERTINIB') %>%
    ggplot(aes(x = Methods, y = `perc_offtarget`, fill = Random)) +
    geom_boxplot(width = 0.2) +
    theme_cowplot() +
    ylim(0, 50) +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = 'none')

plot_small_dispersion <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    filter(Biocontext == 'LNCAP_FORETINIB') %>%
    ggplot(aes(x = Methods, y = `perc_offtarget`, fill = Random)) +
    geom_boxplot(width = 0.2) +
    theme_cowplot() +
    ylim(0, 50) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

right_plots <- ggarrange(plot_bigdispersion, plot_small_dispersion, ncol = 1, nrow = 2) %>% ggpubr::as_ggplot()

disagreement_plots <- ggarrange(disagreement_plot, right_plots, ncol = 2, nrow = 1, widths = c(4, 3), padding = unit(3, 'cm'))



ggsave('plots/panacea_disagreement.svg', plot = disagreement_plots, device = 'svg', width = 13, height = 10, units = 'in', dpi = 200)

disagreement_df %>% group_by(factor(count)) %>%
    summarise(mean = mean(mean_diff, na.rm=TRUE), sd = sd(mean_diff, na.rm = TRUE), 
    min = min(mean_diff, na.rm=TRUE),
    max = max(mean_diff, na.rm=TRUE),
    median = median(mean_diff, na.rm = TRUE))


disagreement_df %>% group_by(factor(count)) %>%
    summarise(mean = mean(sd_diff, na.rm=TRUE), sd = sd(sd_diff, na.rm = TRUE), 
    min = min(sd_diff, na.rm=TRUE),
    max = max(sd_diff, na.rm=TRUE),
    median = median(sd_diff, na.rm = TRUE))

panacea_facet_offtargets <- panacea_data %>% group_by(Biocontext, Methods, Random) %>%
    # create boxplots, real vs random, offtarget_nodes column, facet by bio_context
    ggplot(aes(x = Random, y = `perc_offtarget_nodes`, fill = Random)) +
    geom_boxplot(width = 0.2) +    
    theme_cowplot() +
    theme(strip.text.x = element_text(size = 10)) +
    facet_wrap(~Biocontext)

ggsave('plots/panacea_facet_offtargets.png', plot = panacea_facet_offtargets, device = 'png', width = 13, height = 10, units = 'in', dpi = 200)

panacea_de <- read_tsv('panacea_de_all.tsv')
# per biocontext, count how many genes are DE
panacea_de_count <- panacea_de %>% 
    mutate(bio_context = paste(cell_line, treatment, sep='_')) %>%
    group_by(cell_line, treatment, bio_context) %>%
    filter(padj < 0.05, abs(logFC)> 1.5) %>%
    summarise(count = n()) %>%
    arrange(desc(count)) %>% ungroup()
  
order_celllines <- panacea_de_count %>% group_by(cell_line) %>% summarise(count = sum(count)) %>% arrange(desc(count)) %>% pull(cell_line)
order_treatments <- panacea_de_count %>% group_by(treatment) %>% summarise(count = sum(count)) %>% arrange(count) %>% pull(treatment)

# plot a dotplot of the number of DE genes, cell line in x, treatment in y. The size is the count, the color is if this count is above 30 or below
panacea_de_count_plot <- panacea_de_count %>%
    mutate(cell_line = factor(cell_line, levels = order_celllines),
           treatment = factor(treatment, levels = order_treatments)) %>%
    ggplot(aes(x = cell_line, y = treatment, size = count, color = count > 30)) +
    geom_point() +
    scale_color_manual(values = c('black', 'red')) +
    scale_size_continuous(range = c(1, 10),
                          breaks = c(1, 10, 30, 70, 100, 200, 300), # Define the size breaks you want to show in the legend
                          labels = c("1", "10", "30", "70", "100", "200", "300")) + # Define the labels for the size breaks
    theme_cowplot() +
    # add more size points in legend
    labs(color = '> 30 DE genes') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('plots/panacea_de_threshold.png', plot = panacea_de_count_plot, device = 'png', width = 13, height = 10, units = 'in', dpi = 200)











decryptm_data <- read_csv('decryptm_graphdata_results.csv')

decryptm_data$Random <- ifelse(decryptm_data$Random == 'real', 'False', 'True')

perc_nodes_random <- decryptm_data %>% 
  filter(Random == "True") %>% pull(perc_nodes_with_ec50) %>% mean(., na.rm=TRUE) %>% round(2)

perc_nodes_real <- decryptm_data %>% group_by(Biocontext, Methods, Random) %>%
  filter(Random == "False") %>% pull(perc_nodes_with_ec50) %>% mean(., na.rm=TRUE) %>% round(2)

ec50_plot_random <- decryptm_data %>% group_by(Biocontext, Methods, Random) %>%
  filter(Random == "True") %>%
    # plot histogram of ec50 distribution
    ggplot(aes(x = perc_nodes_with_ec50, pattern=Random)) +
    geom_histogram_pattern(bins = 20, pattern_colour = 'black', fill='grey90', color='black', pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    theme_cowplot() +
    # add a line at mean
    geom_vline(aes(xintercept = mean(perc_nodes_with_ec50, na.rm=TRUE)), linetype = 'dashed') +
    # write annotation with mean and in the 70% of the y axis
    ylim(0, 5000) +
    xlim(0, 40) +    
    annotate(geom = 'text', x = 20, y = 4000, label = paste('Mean:',perc_nodes_random), color = 'black') +
    # min x 0, max undetermined
    theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  axis.text.x = element_blank(),
  axis.title.x = element_blank()
)


ec50_plot_real<- decryptm_data %>% group_by(Biocontext, Methods, Random) %>%
  filter(Random == "False") %>%
    # plot histogram of ec50 distribution
    ggplot(aes(x = perc_nodes_with_ec50, pattern=Random)) +
    geom_histogram_pattern(bins = 20, , pattern_colour = 'black', fill='white', color='black', pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    theme_cowplot() +
    scale_pattern_manual(values = c('False' = 'circle')) +
    # add a line at mean
    geom_vline(aes(xintercept = mean(perc_nodes_with_ec50, na.rm=TRUE)), linetype = 'dashed') +
    # write annotation with mean and in the 70% of the y axis
    ylim(0, 5) +
    annotate(geom = 'text', x = 20, y = 4, label = paste('Mean:',perc_nodes_real), color = 'black') +
    # min x 0, max undetermined
    
    xlim(0, 40) +    
    theme(
  strip.background = element_blank(),
  strip.text.x = element_blank(),
  legend.position = 'right'
)

ec50_plots <- ggarrange(ec50_plot_random, ec50_plot_real, ncol = 1, nrow = 2, heights = c(1, 1))
ggsave('decryptm_ec50.svg', ec50_plots, device = 'svg', width = 13, height = 10, units = 'in', dpi = 200)

ggsave('C:/Users/victo/Desktop/ec50_plot.png', plot = ec50_plot, device = 'png', width = 10, height = 10, units = 'in', dpi = 200)




offtargets <- read_tsv('panacea_offtargets.tsv')
tfs <- read_tsv('tf_activity_results.tsv')
offtargets_list <- offtargets %>% pull(target) %>% unique()
tfs_list <- colnames(tfs) %>% unique()
network_collectri <- read_tsv('network_collectri.sif')
network_collectri_nodes <- union(network_collectri$`# source`, network_collectri$target) %>% unique()
# percentage of offtargets in the network
length(intersect(offtargets_list, network_collectri_nodes)) / length(offtargets_list)
# intersection of the two lists
common <- intersect(offtargets_list, tfs_list)
# create raincloud plot
ggplot(data, aes(x = Methods, y = `offtarget_count`, fill = Random)) +
    geom_boxplot(width = 0.2) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))




ggplot(data, aes(x = Methods, y = `offtarget_count`, fill = Random)) +
    geom_density_ridges(alpha = 0.5) +
    geom_boxplot(width = 0.2) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
    

filtered_data <- decryptm_data %>%
  filter(Random == 'False') %>%
  select(Biocontext, Methods, ec50_values_in_network, ec50_values_out_network)

# Function to convert list-like strings to numeric vectors
convert_to_numeric_vector <- function(list_string) {
  as.numeric(unlist(strsplit(sub("\\[|\\]", "", list_string), ", ")))
}

# Apply the conversion
filtered_data$ec50_values_in_network <- sapply(filtered_data$ec50_values_in_network, convert_to_numeric_vector)
filtered_data$ec50_values_out_network <- sapply(filtered_data$ec50_values_out_network, convert_to_numeric_vector)

# Reshape the data for plotting
ec50_values_long <- filtered_data %>%
  pivot_longer(cols = c('ec50_values_in_network', 'ec50_values_out_network'), names_to = "Location", values_to = "EC50_Value") %>%
  mutate(Location = ifelse(Location == "ec50_values_in_network", "In Network", "Out Network")) %>%
  unnest(EC50_Value)

# numberss
ec50_values_long %>%
  group_by(Biocontext, Methods, Location) %>%
  summarise(mean_ec50 = mean(EC50_Value, na.rm = TRUE), sd_ec50 = sd(EC50_Value, na.rm = TRUE), count = n())%>%
  print(n=100)

ec50_values_long %>%
  group_by(Biocontext) %>%
  summarise(mean_ec50 = mean(EC50_Value, na.rm = TRUE), sd_ec50 = sd(EC50_Value, na.rm = TRUE), count = n())%>%
  print(n=100)

# revisit
ec50_values_long %>%
  group_by(Biocontext, Methods, Location) %>%
  summarise(count = n()) %>%
  ggplot(aes(x = Methods, y = count, fill = Location)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Biocontext) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(ggforce)
library(ggplot2)
library(ggpattern)
ec50_nodes <- ggplot(ec50_values_long) +
  geom_boxplot_pattern(aes(x = Methods, y = EC50_Value, pattern = Location), color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
  labs(x = "Location",
       y = "log EC50 Value") +
theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    facet_wrap(~Biocontext)

ggsave('ec50_nodes.png', ec50_nodes, device = 'png', width = 10, height = 10, units = 'in', dpi = 200)

# numberss

mean_metrics <- panacea_data %>% 
    summarise(`Number of nodes` = mean(`Number of nodes`, na.rm=TRUE), `sd_Number of nodes` = sd(`Number of nodes`, na.rm=TRUE),
    `Number of edges` = mean(`Number of edges`, na.rm=TRUE), `sd_Number of edges` = sd(`Number of edges`, na.rm=TRUE),
    `Connected targets` = mean(`Connected targets`, na.rm=TRUE), `sd_Connected targets` = sd(`Connected targets`, na.rm=TRUE),
    `Perc missing targets` = mean(`Perc missing targets`, na.rm=TRUE), `sd_Perc missing targets` = sd(`Perc missing targets`, na.rm=TRUE),
    `Mean degree centrality` = mean(`Mean degree centrality`, na.rm=TRUE), `sd_Mean degree centrality` = sd(`Mean degree centrality`, na.rm=TRUE),
    `Mean closeness centrality` = mean(`Mean closeness centrality`, na.rm=TRUE), `sd_Mean closeness centrality` = sd(`Mean closeness centrality`, na.rm=TRUE),
    `Mean betweenness centrality` = mean(`Mean betweenness centrality`, na.rm=TRUE), `sd_Mean betweenness centrality` = sd(`Mean betweenness centrality`, na.rm=TRUE))

processed_metrics <- mean_metrics %>% 
  pivot_longer(cols = everything(), names_to = "Metric", values_to = "Mean") %>%
  filter(!is.na(Mean))

panacea_longer_mean <- panacea_data %>% 
    group_by(Methods) %>%
    summarise(`Number of nodes` = mean(`Number of nodes`, na.rm=TRUE), `sd_Number of nodes` = sd(`Number of nodes`, na.rm=TRUE),
    `Number of edges` = mean(`Number of edges`, na.rm=TRUE), `sd_Number of edges` = sd(`Number of edges`, na.rm=TRUE),
    `Connected targets` = mean(`Connected targets`, na.rm=TRUE), `sd_Connected targets` = sd(`Connected targets`, na.rm=TRUE),
    `Mean degree centrality` = mean(`Mean degree centrality`, na.rm=TRUE), `sd_Mean degree centrality` = sd(`Mean degree centrality`, na.rm=TRUE),
    `Mean closeness centrality` = mean(`Mean closeness centrality`, na.rm=TRUE), `sd_Mean closeness centrality` = sd(`Mean closeness centrality`, na.rm=TRUE),
    `Mean betweenness centrality` = mean(`Mean betweenness centrality`, na.rm=TRUE), `sd_Mean betweenness centrality` = sd(`Mean betweenness centrality`, na.rm=TRUE)) %>%
    pivot_longer(cols = -Methods, names_to = "Metric", values_to = "Value") %>%
    filter(!is.na(Value)) %>%
    left_join(processed_metrics, by = c('Metric'))

average_deviation <- panacea_longer_mean %>%
  mutate(ratio = Value / Mean *100 - 100)

library(ComplexHeatmap)

metrics_matrix <- panacea_longer_mean %>% 
  select(Metric, Mean) %>%
  distinct() %>%
  column_to_rownames("Metric") %>%
  as.matrix() %>%
  t()


heatmap_data <- average_deviation %>% select(-Value, -Mean) %>%
  spread(key = Metric, value = ratio)

# Extract the methods and metric names for row and column names
methods <- heatmap_data$Methods
metrics <- colnames(heatmap_data)[!colnames(heatmap_data) %in% "Methods"]

# remove cols in metrics_matrix that are not in metrics, also rearrange the columns to match the order of metrics
metrics_matrix <- metrics_matrix[, metrics] %>% # round values to 2 decimals
  round(2)


# Remove the Methods column to have a pure matrix
heatmap_matrix <- as.matrix(heatmap_data[, metrics])


text_annotation <- HeatmapAnnotation(
    averages = anno_text(
        t(metrics_matrix),  # Transposing the matrix to match the correct orientation
        just = "center", 
        gp = gpar(fontsize = 10, col = "black"),
        border = TRUE,  # Adding a border around the text
        padding = unit(c(2, 2), "mm")  # Adjust padding as needed
    ),
    which = "column",  # Specify that this is a column annotation
    height = unit(2, "cm"),  # Adjust the height of the annotation
    title = "Averages",  # Add a title to the top annotation
    title_gp = gpar(fontsize = 12, fontface = "bold")  # Customize the title appearance
)


text_annotation <- HeatmapAnnotation('Average value across methods' = anno_text(
        t(metrics_matrix), 
        location = 0.5, 
        just = "center",
        rot = 0,
    gp = gpar(col = "black", border = "grey75", fontsize = 11),
    which = "column",
    height = unit(1, "cm"),
    show_name = TRUE
    #change name fontsize
    ),
    annotation_name_gp = gpar(fontsize = 11, fontface = "bold"))

# Draw the heatmap with clustering
heatmap <- Heatmap(heatmap_matrix,
        name = "Deviation from mean (%)",
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        top_annotation = text_annotation,
        rect_gp = gpar(col = "white", lwd = 4),
        show_row_names = TRUE,
        show_column_names = TRUE,
        show_heatmap_legend = FALSE,
        # change axis labels fontsize
        row_names_gp = gpar(fontsize = 11),
        
        column_names_gp = gpar(fontsize = 11),
        col = circlize::colorRamp2(c(-150, 0, 150), c("blue", "grey", "red")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          # 'j' is the column index, 'i' is the row index
          # 'x' and 'y' are the center of the cell
          # 'fill' is the fill color of the cell
          # 'width' and 'height' are the size of the cell
          
          # The value to be placed in the cell
          value <- heatmap_matrix[i, j]
          
          # Formatting the value to a string, you can format it as needed
          text_value <- sprintf("%.2f", value) %>% paste0("%")

          # Adding the text to the cell
          grid.text(text_value, x, y, gp = gpar(fontsize = 10, col = "white", fontface = "bold"))
        })


legend <- Legend(col_fun = circlize::colorRamp2(c(-150, 0, 150), c("blue", "grey", "red")),
                 title = "Deviation from\nmean (%)",
                at = seq(-150, 150, 50),
                # labels = c("-150", "0", "150"),
                 direction = "vertical",
                 title_gp = gpar(fontsize = 11, fontface = "bold"),
                 labels_gp = gpar(fontsize = 10))

# Save the heatmap
png("heatmap.png", width = 9, height = 8, units = "in", res = 300)
draw(heatmap, padding = unit(c(0,0,0,1), 'in'))
draw(legend, x = unit(7.4, "in"), y = unit(1.1, "in"))
dev.off()



mean_metrics %>% View()

panacea_data %>% 
  mutate(Random = ifelse(Random == 'real', 'False', 'True')) %>%
    group_by(Methods, Random) %>%
    summarise(mean = mean(`Number of edges`, na.rm=TRUE), sd = sd(`Number of edges`, na.rm=TRUE)) %>%
    arrange(desc(mean))

plotter_attr <- function(data, attr) {
  p <- data %>% 
    ggplot(aes(x = Methods, y = !!sym(attr), pattern = Random)) +
    geom_boxplot_pattern(color = "black", pattern_fill = "white", pattern_angle = 45, pattern_density = 0.1, pattern_spacing = 0.025, pattern_key_scale_factor = 0.6) +
    # geom_sina(alpha = 0.1, size=0.1) +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Location",
         y = attr) +
         facet_wrap(~`Study ID`)

  ggsave(paste0('C:/Users/victo/Desktop/boxplot_', attr, '.png'), plot = p, device = 'png', width = 10, height = 10, units = 'in', dpi = 200)
}

panacea_data <- panacea_data %>% rename("perc_offtargets" = "perc_offtarget",
                        "perc_offtargets_nodes" = "perc_offtarget_nodes",
                        "perc_offtargets_edges" = "perc_offtarget_edges")

full_data <- bind_rows(panacea_data, decryptm_data)

panacea_data %>% colnames()
plotter_attr(full_data, 'Number of nodes')
plotter_attr(full_data, 'Number of edges')
plotter_attr(full_data, 'Connected targets')
plotter_attr(full_data, 'Perc missing targets')
plotter_attr(full_data, 'Mean degree centrality')
plotter_attr(full_data, 'Mean closeness centrality')
plotter_attr(full_data, 'Mean betweenness centrality')
plotter_attr(full_data, 'perc_offtargets_nodes')
plotter_attr(full_data, 'perc_offtargets')
