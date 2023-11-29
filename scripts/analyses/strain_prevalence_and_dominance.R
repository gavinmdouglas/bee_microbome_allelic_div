rm(list = ls(all.names = TRUE))

# Investigation into where specific strains consistently dominant/minor?
# Summarize first what the relative abundance of the most abundant strain is across samples/species
# Then investigate whether certain strains are consistently dominant minor

# Do this for both StrainFacts and StrainGST-identified strains.

strainfacts_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds')
straingst_relabun <- readRDS('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/strainge/RDS/straingst_relabun.rds')

# Make sure the same samples are compared.
strain_species <- intersect(names(strainfacts_relabun), names(straingst_relabun))
strain_samples <- list()

for (sp in strain_species) {
  tmp_datasets <- intersect(names(strainfacts_relabun[[sp]]), names(straingst_relabun[[sp]]))
  if (length(tmp_datasets) == 0) { next }
  
  for (tmp_d in tmp_datasets) {
    tmp_intersecting_samp <- intersect(colnames(strainfacts_relabun[[sp]][[tmp_d]]), colnames(straingst_relabun[[sp]][[tmp_d]]))
    if (length(tmp_intersecting_samp) == 0) { next }
    
    if (! sp %in% names(strain_samples))   {
      strain_samples[[sp]] <- list()
    }
    
    strain_samples[[sp]][[tmp_d]] <- tmp_intersecting_samp
  }
}

# Get prevalence per strain across samples,
# mean relative abundance, when present,
# and proportion of times it is the most frequent strain.

raw_out_strainfacts <- list()
raw_out_straingst <- list()

# Also run multinomial tests on
# 1) Strains being equally likely to occur across samples.
# 2) Strains being equally likely to be the most common (when they are present).
raw_out <- list()

for (sp in names(strain_samples)) {
  
  for (d in names(strain_samples[[sp]])) {

    strainfacts_tab <- strainfacts_relabun[[sp]][[d]]
    strainfacts_tab_no_zero <- strainfacts_tab
    strainfacts_tab_no_zero[strainfacts_tab_no_zero == 0] <- NA

    strainfacts_summary <- data.frame(Species = sp,
                                      Dataset = d,
                                      Tool = 'StrainFacts',
                                      Strain = rownames(strainfacts_tab),
                                      prevalence = rowSums(strainfacts_tab > 0),
                                      prevalence_prop = rowSums(strainfacts_tab > 0) / ncol(strainfacts_tab),
                                      prevalence_norm = rowSums(strainfacts_tab > 0) / sum(rowSums(strainfacts_tab > 0)),
                                      prevalence_norm_over_exp = NA,
                                      mean_relabun_present = rowMeans(strainfacts_tab_no_zero, na.rm = TRUE),
                                      times_most_frequent = 0)
    rownames(strainfacts_summary) <- strainfacts_summary$Strain

    strain_most_frequent <- table(sapply(strainfacts_tab, function(x) { rownames(strainfacts_tab)[which.max(x)] }))

    strainfacts_summary[names(strain_most_frequent), 'times_most_frequent'] <- strain_most_frequent
    strainfacts_summary$prop_most_frequent <- strainfacts_summary$times_most_frequent / ncol(strainfacts_tab)
    strainfacts_summary$prop_present_most_frequent <- strainfacts_summary$times_most_frequent / strainfacts_summary$prevalence

    # Test for whether some strains are more likely to be present than others.
    if (nrow(strainfacts_tab) > 1) { 
      strainfacts_exp_prev_prob_vec <- rep(1 / nrow(strainfacts_tab), nrow(strainfacts_tab))
      strainfacts_multinomial_prev_p <- XNomial::xmonte(strainfacts_summary$prevalence_norm,
                                                        strainfacts_exp_prev_prob_vec,
                                                        detail = 0)$pLLR

      strainfacts_max_prev_ratio <- max(strainfacts_summary$prevalence_norm) / strainfacts_exp_prev_prob_vec[1]
      strainfacts_min_prev_ratio <- min(strainfacts_summary$prevalence_norm) / strainfacts_exp_prev_prob_vec[1]
      strainfacts_summary$prevalence_norm_over_exp <- strainfacts_summary$prevalence_norm / strainfacts_exp_prev_prob_vec[1]
      # Test whether some strains are more likely to be dominant (when present).
      strainfacts_multinomial_dom_p <- XNomial::xmonte(strainfacts_summary$times_most_frequent,
                                                       strainfacts_summary$prevalence_prop,
                                                       detail = 0)$pLLR
      
      strainfacts_dom_ratio <- ((strainfacts_summary$times_most_frequent + 1) /
                                  ((sum(strainfacts_summary$times_most_frequent)) + length(strainfacts_summary$times_most_frequent))) /
                   (strainfacts_summary$prevalence_prop / sum(strainfacts_summary$prevalence_prop))

      raw_out[[paste(sp, d, 'strainfacts')]] <- data.frame(Species = sp,
                                                           Dataset = d,
                                                           Tool = 'StrainFacts',
                                                           exp_prop_prev = strainfacts_exp_prev_prob_vec[1],
                                                           num_strains = nrow(strainfacts_tab),
                                                           num_samples = ncol(strainfacts_tab),
                                                           mean_strains_per_sample = mean(colSums(strainfacts_tab > 0)),
                                                           max_strains_per_sample = max(colSums(strainfacts_tab > 0)),

                                                           mean_prev = mean(strainfacts_summary$prevalence_prop),
                                                           prev_multinomial_p = strainfacts_multinomial_prev_p,
                                                           max_prev_ratio = strainfacts_max_prev_ratio,
                                                           min_prev_ratio = strainfacts_min_prev_ratio,
                                                           
                                                           multinomial_dom_p = strainfacts_multinomial_dom_p,
                                                           max_dom_ratio = strainfacts_dom_ratio[which.max(strainfacts_dom_ratio)],
                                                           min_dom_ratio = strainfacts_dom_ratio[which.min(strainfacts_dom_ratio)])
      
    }

    straingst_tab <- straingst_relabun[[sp]][[d]]
    straingst_tab_no_zero <- straingst_tab
    straingst_tab_no_zero[straingst_tab_no_zero == 0] <- NA
    
    straingst_summary <- data.frame(Species = sp,
                                      Dataset = d,
                                      Tool = 'StrainGST',
                                      Strain = rownames(straingst_tab),
                                      prevalence = rowSums(straingst_tab > 0),
                                      prevalence_prop = rowSums(straingst_tab > 0) / ncol(straingst_tab),
                                      prevalence_norm = rowSums(straingst_tab > 0) / sum(rowSums(straingst_tab > 0)),
                                      prevalence_norm_over_exp = NA,
                                      mean_relabun_present = rowMeans(straingst_tab_no_zero, na.rm = TRUE),
                                      times_most_frequent = 0)
    rownames(straingst_summary) <- straingst_summary$Strain
    
    strain_most_frequent <- table(sapply(straingst_tab, function(x) { rownames(straingst_tab)[which.max(x)] }))
    
    straingst_summary[names(strain_most_frequent), 'times_most_frequent'] <- strain_most_frequent
    straingst_summary$prop_most_frequent <- straingst_summary$times_most_frequent / ncol(straingst_tab)
    straingst_summary$prop_present_most_frequent <- straingst_summary$times_most_frequent / straingst_summary$prevalence
    
    if (nrow(straingst_tab) > 1) {
      # Test for whether some strains are more likely to be present than others.
      straingst_exp_prev_prob_vec <- rep(1 / nrow(straingst_tab), nrow(straingst_tab))
      straingst_multinomial_prev_p <- XNomial::xmonte(straingst_summary$prevalence,
                                                        straingst_exp_prev_prob_vec,
                                                        detail = 0)$pLLR

      straingst_max_prev_ratio <- max(straingst_summary$prevalence_norm) / straingst_exp_prev_prob_vec[1]
      straingst_min_prev_ratio <- min(straingst_summary$prevalence_norm) / straingst_exp_prev_prob_vec[1]
      straingst_summary$prevalence_norm_over_exp <- straingst_summary$prevalence_norm / straingst_exp_prev_prob_vec[1]
    
      # Test whether some strains are more likely to be dominant (when present).
      straingst_multinomial_dom_p <- XNomial::xmonte(straingst_summary$times_most_frequent,
                                                       straingst_summary$prevalence_prop,
                                                       detail = 0)$pLLR
      
      straingst_dom_ratio <- ((straingst_summary$times_most_frequent + 1) /
                                ((sum(straingst_summary$times_most_frequent)) + length(straingst_summary$times_most_frequent))) /
                             (straingst_summary$prevalence_prop / sum(straingst_summary$prevalence_prop))
      
      raw_out[[paste(sp, d, 'straingst')]] <- data.frame(Species = sp,
                                                           Dataset = d,
                                                           Tool = 'straingst',
                                                           exp_prop_prev = straingst_exp_prev_prob_vec[1],
                                                           num_strains = nrow(straingst_tab),
                                                           num_samples = ncol(straingst_tab),
                                                           mean_strains_per_sample = mean(colSums(straingst_tab > 0)),
                                                           max_strains_per_sample = max(colSums(straingst_tab > 0)),
                                                           
                                                           mean_prev = mean(straingst_summary$prevalence_prop),
                                                           prev_multinomial_p = straingst_multinomial_prev_p,
                                                           max_prev_ratio = straingst_max_prev_ratio,
                                                           min_prev_ratio = straingst_min_prev_ratio,
                                                           
                                                           multinomial_dom_p = straingst_multinomial_dom_p,
                                                           max_dom_ratio = straingst_dom_ratio[which.max(straingst_dom_ratio)],
                                                           min_dom_ratio = straingst_dom_ratio[which.min(straingst_dom_ratio)])
    }
    
    raw_out_strainfacts[[paste(sp, d, 'strainfacts')]] <- strainfacts_summary
    raw_out_straingst[[paste(sp, d, 'straingst')]] <- straingst_summary

  }
}

combined_multinomial_out <- do.call(rbind, raw_out)

combined_multinomial_out$Sig_prev <- 'No'
combined_multinomial_out[which(combined_multinomial_out$prev_multinomial_p < 0.05), 'Sig_prev'] <- 'Yes'

combined_multinomial_out$Sig_dom <- 'No'
combined_multinomial_out[which(combined_multinomial_out$multinomial_dom_p < 0.05), 'Sig_dom'] <- 'Yes'

strainfacts_multinomial_out <- combined_multinomial_out[which(combined_multinomial_out$Tool == 'StrainFacts'), ]
straingst_multinomial_out <- combined_multinomial_out[which(combined_multinomial_out$Tool == 'straingst'), ]

ggplot(data = strainfacts_multinomial_out,
       aes(x = log2(max_prev_ratio),
           y = log2(min_prev_ratio),
           colour = Sig_prev)) +
  geom_point(aes(pch=Dataset), size = 4) +
  facet_wrap(Species ~ .) +
  xlim(c(0, 6)) +
  ylim(c(-6, 0)) +
  geom_abline(slope = -1,
              intercept = 0,
              color='black',
              lty=2,
              lwd=1) +
  xlab('log2(Max. prevalence ratio)') +
  ylab('log2(Min. prevalence ratio)') +
  theme_bw() +
  scale_colour_manual(values=c('grey75', 'cornflowerblue'))
  

strainfacts_multinomial_out_dominance <- strainfacts_multinomial_out[which(strainfacts_multinomial_out$max_strains_per_sample > 1), ]

ggplot(data = strainfacts_multinomial_out_dominance,
       aes(x = log2(max_dom_ratio),
           y = log2(min_dom_ratio),
           colour = Sig_dom)) +
  geom_point(aes(pch=Dataset), size = 4) +
  facet_wrap(Species ~ .) +
  xlim(c(0, 6)) +
  ylim(c(-6, 0)) +
  geom_abline(slope = -1,
              intercept = 0,
              color='black',
              lty=2,
              lwd=1) +
  xlab('log2(Max. dominance ratio)') +
  ylab('log2(Min. dominance ratio)') +
  theme_bw() +
  scale_colour_manual(values=c('grey75', 'cornflowerblue'))



straingst_multinomial_out_dominance <- straingst_multinomial_out[which(straingst_multinomial_out$max_strains_per_sample > 1), ]

ggplot(data = straingst_multinomial_out_dominance,
       aes(x = log2(max_dom_ratio),
           y = log2(min_dom_ratio),
           colour = Sig_dom)) +
  geom_point(aes(pch=Dataset), size = 4) +
  facet_wrap(Species ~ .) +
  xlim(c(0, 6)) +
  ylim(c(-6, 0)) +
  geom_abline(slope = -1,
              intercept = 0,
              color='black',
              lty=2,
              lwd=1) +
  xlab('log2(Max. dominance ratio)') +
  ylab('log2(Min. dominance ratio)') +
  theme_bw() +
  scale_colour_manual(values=c('grey75', 'cornflowerblue'))

ggplot(data = straingst_multinomial_out,
       aes(x = log2(max_prev_ratio),
           y = log2(min_prev_ratio),
           colour = Sig_prev)) +
  geom_point(aes(pch=Dataset), size = 4) +
  facet_wrap(Species ~ .) +
  xlim(c(0, 6)) +
  ylim(c(-6, 0)) +
  geom_abline(slope = -1,
              intercept = 0,
              color='black',
              lty=2,
              lwd=1) +
  xlab('log2(Max. prevalence ratio)') +
  ylab('log2(Min. prevalence ratio)') +
  theme_bw() +
  scale_colour_manual(values=c('grey75', 'cornflowerblue'))

full_strainfacts_summary <- do.call(rbind, raw_out_strainfacts)
rownames(full_strainfacts_summary) <- NULL

full_strainfacts_summary_atleast5 <- full_strainfacts_summary[which(full_strainfacts_summary$prevalence >= 5), ]

ggplot(data = full_strainfacts_summary,
       aes(x = prevalence_prop, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Proportion of samples with strain')

ggplot(data = full_strainfacts_summary,
       aes(x = prevalence_norm, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Proportion of times strain occurred out of total occurrences')

ggplot(data = full_strainfacts_summary,
       aes(x = prevalence_norm_over_exp, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Prevalence ratio: (Proportion of times strain occurred out of total occurrences) / (Expected proportion under equality)')


ggplot(data = full_strainfacts_summary,
       aes(x = prop_most_frequent, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Proportion of times strain is dominant (whether present or not)')

ggplot(data = full_strainfacts_summary,
       aes(x = prop_present_most_frequent, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Proportion of times strain is dominant (when present)')

ggplot(data = full_strainfacts_summary,
       aes(x = prevalence_norm, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Proportion of times strain occurred out of total occurrences')

ggplot(data = full_strainfacts_summary,
       aes(x = prevalence_norm_over_exp, y = Dataset)) +
  geom_violin() +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(Species ~ .) +
  xlab('Prevalence ratio: (Proportion of times strain occurred out of total occurrences) / (Expected proportion under equality)')


