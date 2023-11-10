rm(list = ls(all.names = TRUE))

library(bestNormalize)
library(glmmTMB)
library(performance)

which.median <- function(x) which.min(abs(x - median(x)))

breakdown <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/statistics/num_strain_vs_num_allele_breakdown.tsv',
                        header = TRUE, sep = '\t', stringsAsFactors = FALSE)
breakdown <- breakdown[-which(is.na(breakdown$cor_present_p)), ]
breakdown <- breakdown[which(breakdown$sample_prev >= 30), ]
breakdown <- breakdown[which(breakdown$num_alleles >= 10), ]
breakdown <- breakdown[which(breakdown$Species != 'Commensalibacter_sp'), ]

breakdown$tau_present_norm <- bestNormalize::orderNorm(breakdown$tau_present)$x

species_mean_val <- aggregate(x = tau_present ~ Species, data = breakdown, FUN = mean)
median_species <- species_mean_val$Species[which.median(species_mean_val$tau_present)]

breakdown$Species <- factor(breakdown$Species,
                            levels = c(median_species, species_mean_val$Species[-which(species_mean_val$Species == median_species)]))

lm_out <- lm(tau_present ~ Species, data = breakdown)

model_out <- glmmTMB::glmmTMB(formula = tau_present ~ Species + (1 | num_alleles) + (1 | sample_prev),
                               data = breakdown,
                               se = TRUE,
                               family = 'gaussian',
                               control = glmmTMBControl(optimizer = nlminb,
                                                        parallel = 40,
                                                        profile = TRUE,
                                                        optCtrl = list(iter.max = 1000,
                                                                       eval.max = 1000)))

model_out_summary <- summary(model_out)
species_coefficients <- data.frame(model_out_summary$coefficients$cond[grep('^Species', rownames(model_out_summary$coefficients$cond)), ])
species_coefficients$variable <- gsub('^Species', '', rownames(species_coefficients))
species_coefficients$variable <- gsub('_', ' ', species_coefficients$variable)
species_coefficients$variable <- gsub(' sp$', ' sp.', species_coefficients$variable)
species_coefficients$variable <- gsub('Bifidobacterium coryneforme indicum', 'Bifidobacterium cor./indicum', species_coefficients$variable)

species_coefficients$P_category <- 'P >= 0.05'
species_coefficients$P_category[which(species_coefficients$Pr...z.. < 0.05)] <- 'P < 0.05'

species_coefficients$variable <- factor(species_coefficients$variable,
                                        levels = species_coefficients$variable[order(species_coefficients$Estimate)])

saveRDS(object = species_coefficients, file = '/home/gdouglas/tmp/species_coefficients.rds')

species_coefficients_barplot <- ggplot(data = species_coefficients, aes(x = Estimate, y = variable, fill = P_category)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = c('firebrick', 'grey85')) +
  theme_bw() +
  ylab('Species') +
  xlab('Intercept coefficient') +
  theme(axis.text.y = element_text(face = 'italic')) +
  geom_vline(xintercept = 0, linetype="dotted", 
             color = "black") +
  geom_errorbar(aes(xmin = Estimate - Std..Error,
                    xmax = Estimate + Std..Error),
                width = 0.2, color = "black") +
  guides(fill = guide_legend('Coefficient\ntest result'))


