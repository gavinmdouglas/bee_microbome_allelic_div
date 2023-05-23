rm(list = ls(all.names = TRUE))

library(ape)
library(ggplot2)
library(ggtree)
library(cowplot)

# Apilactobacillus apinoruim vs kunkeei
Apilactobacillus_phylo <- read.tree(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Apilactobacillus.tree')
Apilactobacillus_phylo_plot <- ggplot(Apilactobacillus_phylo, aes(x, y)) +
                                      geom_tree() +
                                      theme_tree2() +
                                      geom_tiplab(as_ylab = TRUE, color = 'black') +
                                      geom_rootedge() +
                                      ggtitle(expression(italic("Apilactobacillus"))) +
                                      theme(plot.title = element_text(hjust = 0.5))

# Bartonella_apis
Bartonella_apis_phylo <- read.tree(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Bartonella_apis.tree')
Bartonella_apis_phylo_plot <- ggplot(Bartonella_apis_phylo, aes(x, y)) +
                                      geom_tree() +
                                      theme_tree2() +
                                      geom_tiplab(as_ylab = TRUE, color = 'black') +
                                      geom_rootedge() +
                                      ggtitle(expression(italic("Bartonella apis"))) +
                                      theme(plot.title = element_text(hjust = 0.5))

# All Bifidobacterium
bifido_phylo <- read.tree(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Bifidobacterium.tree')
bifido_phylo_plot <- ggplot(bifido_phylo, aes(x, y)) +
                            geom_tree() +
                            theme_tree2() +
                            geom_tiplab(as_ylab = TRUE, color = 'black') +
                            geom_rootedge() +
                            ggtitle(expression(italic("Bifidobacterium"))) +
                            theme(plot.title = element_text(hjust = 0.5))



# Bombella apis / Bombella sp / Parasaccharibacter_apium / Saccharibacter_sp
bombella_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Bombella_and_related.tree')
bombella_phylo_plot <- ggplot(bombella_phylo, aes(x, y)) +
                              geom_tree() +
                              theme_tree2() +
                              geom_tiplab(as_ylab = TRUE, color = 'black') +
                              geom_rootedge() +
                              ggtitle(expression(italic("Bombella"))) +
                              theme(plot.title = element_text(hjust = 0.5))

# Bombilactobacillus mellis / mellifer
bombilactobacillus_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Bombilactobacillus_and_related.tree')

bombilactobacillus_phylo_plot <- ggplot(bombilactobacillus_phylo, aes(x, y)) +
                                    geom_tree() +
                                    theme_tree2() +
                                    geom_tiplab(as_ylab = TRUE, color = 'black') +
                                    geom_rootedge() +
                                    ggtitle(expression(italic("Bombilactobacillus"))) +
                                    theme(plot.title = element_text(hjust = 0.5))

# Commensalibacter_sp
Commensalibacter_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Commensalibacter_sp.tree')
Commensalibacter_phylo_plot <- ggplot(Commensalibacter_phylo, aes(x, y)) +
                                      geom_tree() +
                                      theme_tree2() +
                                      geom_tiplab(as_ylab = TRUE, color = 'black') +
                                      geom_rootedge() +
                                      ggtitle(expression(italic("Commensalibacter sp."))) +
                                      theme(plot.title = element_text(hjust = 0.5))


# Frischella_perrara / Gilliamella apis / Gilliamella apicola
gamma_proteo_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Gilliamella_and_Frishchella_perrara.tree')

gamma_proteo_phylo_plot <- ggplot(gamma_proteo_phylo, aes(x, y)) +
                                  geom_tree() +
                                  theme_tree2() +
                                  geom_tiplab(as_ylab = TRUE, color = 'black') +
                                  geom_rootedge() +
                                  ggtitle('Gammaproteobacteria') +
                                  theme(plot.title = element_text(hjust = 0.5))



# Lactobacilli
lactobacilli_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Lactobacillus.tree')

Lactobacillus_sp_to_exclude <- c("GCA_000760615.1", "GCF_016101625.1", "GCA_003693025.1",
                                 "GCF_016100885.1", "GCF_016100935.1", "GCF_016100975.1")

tip_label_colours <- rep(x = 'black', times = length(lactobacilli_phylo$tip.label))
tip_label_colours[which(gsub('^.*\\|', '', lactobacilli_phylo$tip.label) %in% Lactobacillus_sp_to_exclude)] <- 'red'
  
lactobacilli_phylo_plot <- ggplot(lactobacilli_phylo, aes(x, y)) +
                                  geom_tree() +
                                  theme_tree2() +
                                  geom_tiplab(as_ylab = TRUE, colour = tip_label_colours) +
                                  geom_rootedge() +
                                  ggtitle(expression(italic("Lactobacillus"))) +
                                  theme(plot.title = element_text(hjust = 0.5))

# Serratia_marcescens
Serratia_marcescens_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Serratia_marcescens.tree')
Serratia_marcescens_phylo_plot <- ggplot(Serratia_marcescens_phylo, aes(x, y)) +
                                    geom_tree() +
                                    theme_tree2() +
                                    geom_tiplab(as_ylab = TRUE) +
                                    geom_rootedge() +
                                    ggtitle(expression(italic("Serratia marcescens"))) +
                                    theme(plot.title = element_text(hjust = 0.5))


# Snodgrassella_alvi
Snodgrassella_alvi_phylo <- read.tree('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/ANI_NJ_trees/Snodgrassella_alvi.tree')
Snodgrassella_alvi_phylo_plot <- ggplot(Snodgrassella_alvi_phylo, aes(x, y)) +
  geom_tree() +
  theme_tree2() +
  geom_tiplab(as_ylab = TRUE) +
  geom_rootedge() +
  ggtitle(expression(italic("Snodgrassella alvi"))) +
  theme(plot.title = element_text(hjust = 0.5))

multi.species_combined_ANI_NJ_plot <- plot_grid(lactobacilli_phylo_plot, bifido_phylo_plot, gamma_proteo_phylo_plot,
                                           Apilactobacillus_phylo_plot, bombella_phylo_plot, bombilactobacillus_phylo_plot,
                                           rel_heights = c(1, 0.5),
                                           labels = letters[1:6])

single.species_combined_ANI_NJ_plot <- plot_grid(Bartonella_apis_phylo_plot,
                                             Commensalibacter_phylo_plot,
                                             Serratia_marcescens_phylo_plot,
                                             Snodgrassella_alvi_phylo_plot,
                                             labels = letters[1:4])



ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_ANI_NJ_multi.species_trees.pdf',
       plot = multi.species_combined_ANI_NJ_plot,
       device = 'pdf',
       dpi = 400,
       width = 20,
       height = 12)

ggsave(filename = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Supp_ANI_NJ_single.species_trees.pdf',
       plot = single.species_combined_ANI_NJ_plot,
       device = 'pdf',
       dpi = 400,
       width = 10,
       height = 9)
