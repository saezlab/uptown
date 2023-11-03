library(ggplot2)
library(ggridges)
library(tidyverse)

# create example data
data <- read_csv('panacea_graphdata_results.csv')
data$Random %>% unique() %>% length()
# change the random column so that it contains True when Random==real and False otherwise
data$Random <- ifelse(data$Random == 'real', 'False', 'True')

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
    




