rm(list = ls(all.names = TRUE))

library('MASS')
library('reshape2')

strain_abun <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/core/RDS/strainfacts_core.genome_comm.rds")

all_genes_abun <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds")

mean_depth <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_mean_depth.tsv.gz',
                         header = TRUE, row.names = 1, sep = '\t', stringsAsFactors = FALSE)

species <- 'Bifidobacterium_asteroides'
dataset <- 'Ellegaard2020'
gene <- 'Bifidobacterium_asteroides_eda'
strain_tab <- strain_abun[[species]][[dataset]]

allele_tab <- all_genes_abun[[dataset]][[species]][[gene]]
allele_tab$strain <- paste('allele_', allele_tab$strain, sep = '')
allele_tab <- reshape2::dcast(data = allele_tab, formula = sample ~ strain, fill = 0, value.var = 'community')
rownames(allele_tab) <- allele_tab$sample
allele_tab <- allele_tab[, -which(colnames(allele_tab) == 'sample'), drop = FALSE]

intersecting_samples <- intersect(colnames(strain_tab), rownames(allele_tab))
strain_tab <- strain_tab[, intersecting_samples, drop = FALSE]
allele_tab <- allele_tab[intersecting_samples, , drop = FALSE]

strain_tab[strain_tab > 0] <- 1
allele_tab[allele_tab > 0] <- 1

strain_tab_pinv <- ginv(t(as.matrix(strain_tab)))
X <- strain_tab_pinv %*% as.matrix(allele_tab)

t(as.matrix(strain_tab)) %*% X 

X_clean <- X
X_clean[X_clean < 0] <- 0
for (i in 1:nrow(X_clean)) {
  
}


# Mean depth of variant sites
mean_variant_depth <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/pandora_running/pandora_out_illumina/pandora_mean_variant_site_depth.tsv",
                                 header = TRUE, sep = "\t", row.names = 1)

# Test genes with tractable #s of combos to test
# 44 Snodgrassella_group_660 11
# 60 Snodgrassella_group_629  7
# 68        Gilliamella_sixA  7
# 75  Gilliamella_group_1302  7
# 80     Firm5_pepO_2~~~pepO 30
# 81        Firm5_group_3578 18
# 83      Snodgrassella_fkpA 15
# 88        Firm5_group_1574 12

gene <- "Snodgrassella_fkpA"
num_orig_haplotypes <- 15

gene_samples_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/prepped_input/", 
                           gene, "_samples.txt", sep = "")
gene_samples <- read.table(file = gene_samples_path, stringsAsFactors = FALSE)$V1

phylotype <- strsplit(gene, split = "_")[[1]][1]
phylotype_strains <- core_haplotype_abun[[phylotype]]

otu_tab_path <- paste("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/StrainFinder_running/test_representative_genes/output/pandora_prepped_w_struct/otu_tables/otu_table.",
                      gene, ".", as.character(num_orig_haplotypes), ".txt", sep = "")
otu_tab <- read.table(file = otu_tab_path, header = FALSE, skip = 1, sep = "\t")
otu_tab[otu_tab < 0.01] <- 0

rownames(otu_tab) <- gene_samples

# Remove haplotypes that never intersect with any strains.
gene_and_strain_samples <- gene_samples[which(gene_samples %in% rownames(phylotype_strains))]
intersecting_haplotypes_i <- which(colSums(otu_tab[gene_and_strain_samples, ]) > 0)
otu_tab <- otu_tab[, intersecting_haplotypes_i, drop = FALSE]

num_haplotypes <- ncol(otu_tab)

# Add missing samples to the strain and OTU tables (with only abundances given by the pseudocounts below).
all_samples_to_test <- unique(c(gene_samples, rownames(phylotype_strains)))

gene_samples_with_no_strain <- gene_samples[which(! gene_samples %in% rownames(phylotype_strains))]
if (length(gene_samples_with_no_strain) > 0) {
  phylotype_strains[gene_samples_with_no_strain, ] <- 0
}

strain_samples_with_no_gene <- rownames(phylotype_strains)[which(! rownames(phylotype_strains) %in% gene_samples)]
if (length(strain_samples_with_no_gene) > 0) {
  otu_tab[strain_samples_with_no_gene, ] <- 0
}

# For each strain, figure out which haplotypes it could ostensibly encode
# The clearest cases to exclude are those where the strain and haplotype
# never co-occur in the same sample. Can also ignore strain / haplotype
# possibilities that don't co-occur at least 75% of the time when the strain is there.

possible_strain_haplotypes <- list()

strains2ignore <- c()

for (strain_i in 1:ncol(phylotype_strains)) {
  
  haplotypes_to_add <- c(1:ncol(otu_tab), ncol(otu_tab) + 1)
  
  samples_w_strain <- rownames(phylotype_strains)[which(phylotype_strains[, strain_i] > 0)]
  
  sample_haplotype_to_ignore <- which(colSums(otu_tab[samples_w_strain, , drop = FALSE] > 0) / length(samples_w_strain) < 0.75)
  
  if (length(sample_haplotype_to_ignore) > 0) {
    if (length(haplotypes_to_add[-sample_haplotype_to_ignore]) == 1 && haplotypes_to_add[-sample_haplotype_to_ignore] == ncol(otu_tab) + 1) {
      strains2ignore <- c(strains2ignore, strain_i)
    }
    
    possible_strain_haplotypes[[strain_i]] <- haplotypes_to_add[-sample_haplotype_to_ignore]
  } else {
    possible_strain_haplotypes[[strain_i]] <- haplotypes_to_add
  }
}

if (length(strains2ignore) > 0) {
  possible_strain_haplotypes <- possible_strain_haplotypes[-strains2ignore]
  phylotype_strains <- phylotype_strains[, -strains2ignore]
}

remaining_haplotypes <- c()
for (i in 1:length(possible_strain_haplotypes)) {
  remaining_haplotypes <- c(remaining_haplotypes, possible_strain_haplotypes[[i]])
}
remaining_haplotypes <- unique(remaining_haplotypes)

possible_strain_haplotype_combos <- expand.grid(possible_strain_haplotypes)

# possible_strain_haplotype_combos_all_represented <- possible_strain_haplotype_combos
# 
# 
# for (i in remaining_haplotypes) {
#   if (i == num_haplotypes + 1) { next }
#   possible_strain_haplotype_combos_all_represented <- possible_strain_haplotype_combos_all_represented[-which(rowSums(possible_strain_haplotype_combos_all_represented == i) == 0), ]
# }
# 
# # Remove all combos where a strain is assigned a haplotype it never intersect with.
# for (strain_i in 1:ncol(phylotype_strains)) {
# 
#   samples_w_strain <- rownames(phylotype_strains)[which(phylotype_strains[, strain_i] > 0)]
#   
#   samples_w_strain_w_genes <- samples_w_strain[which(samples_w_strain %in% gene_samples)]
#   
#   sample_haplotype_missing <- which(sapply(otu_tab[samples_w_strain_w_genes, , drop = FALSE], max) == 0)
#   
#   if (length(sample_haplotype_missing) > 0) {
#     for (strain_haplotype_i in sample_haplotype_missing) {
#       rows2rm <- which(possible_strain_haplotype_combos[, strain_i] == strain_haplotype_i)
#       
#       if (length(rows2rm) > 0) {
#         possible_strain_haplotype_combos <- possible_strain_haplotype_combos[-rows2rm, ]
#       }
#     }
#   }
# }

# Add pseudocount of 0.01 to haplotype and strain abundances for all 0s.
phylotype_strains[phylotype_strains == 0] <- 0.01
otu_tab[otu_tab == 0] <- 0.01

non_gene_samples <- rownames(phylotype_strains)[which(! rownames(phylotype_strains) %in% gene_samples)]

# Compute depth values to use in advance.
gene_mean_variant_depth <- list()
for (gene_sample in gene_samples) {
  gene_mean_variant_depth[[gene_sample]] <- round(mean_variant_depth[gene, gene_sample])
  if (gene_mean_variant_depth[[gene_sample]] < 10) {
    stop("For some reason mean depth is lower than 10 reads...")
  }
}

other_mean_variant_depth <- list()
mean_variant_depth_phylotype <- mean_variant_depth[grep(phylotype, rownames(mean_variant_depth)), ]
for (non_gene_sample in all_samples_to_test[which(! all_samples_to_test %in% gene_samples)]) {
  present_phylotype_genes <- rownames(mean_variant_depth_phylotype)[which(mean_variant_depth_phylotype[, non_gene_sample] > 0)]
  other_mean_variant_depth[[non_gene_sample]] <- round(mean(mean_variant_depth[present_phylotype_genes, non_gene_sample]))
  if (other_mean_variant_depth[[non_gene_sample]] < 10) {
    stop("For some reason mean depth is lower than 10 reads...")
  }
}


compute_haplotype_exp_prob_based_on_strains <- function(haplotype_dist, strain_abun, num_hap, min_prob = 0.01) {
  
  haplotype_prob <- rep(min_prob, num_hap)
  
  for(strain_present in unique(haplotype_dist)) {
    
    # Skip "haplotype" that indicates a state of not encoding haplotype.
    if (strain_present == num_hap + 1) { next }
    
    haplotype_prob[strain_present] <- sum(strain_abun[which(haplotype_dist == strain_present)])
  }
  
  return(haplotype_prob)
}


# Get observed haplotype abundances calculated in advance.

sample_obs_abun <- list()
for (gene_sample in gene_samples) {
  sample_obs_abun[[gene_sample]] <- round(gene_mean_variant_depth[[gene_sample]] * as.numeric(otu_tab[gene_sample, ]))
  
  # if (length(which(sample_obs_abun[[gene_sample]] == 0)) > 0) {
  #   sample_obs_abun[[gene_sample]][which(sample_obs_abun[[gene_sample]] == 0)] <- 1
  # }
}

for (non_gene_sample in non_gene_samples) {
  sample_obs_abun[[non_gene_sample]] <- round(other_mean_variant_depth[[non_gene_sample]] * rep(1, ncol(otu_tab)))
}


combo_loglik <- function(param) {
  
  # Discretize parameter values.
  possible_strain_haplotype_combos <- floor(param[1:ncol(phylotype_strains)])
  
  #possible_strain_haplotype_combos <- possible_strain_haplotype_combos[1, ]
  
  # Make sure that min param value is 1 and the max value is the num haplotypes + 1
  if (length(which(possible_strain_haplotype_combos == 0)) > 0) {
    possible_strain_haplotype_combos[which(possible_strain_haplotype_combos == 0)] <- 1
  }
  
  if (length(which(possible_strain_haplotype_combos > (num_haplotypes + 1))) > 0) {
    possible_strain_haplotype_combos[which(possible_strain_haplotype_combos > (num_haplotypes + 1))] <- num_haplotypes + 1
  }
  
  combo_likelihood <- 0
  
  for (gene_sample in gene_samples) {
    
    exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = as.numeric(possible_strain_haplotype_combos),
                                                            strain_abun = as.numeric(phylotype_strains[gene_sample, ]),
                                                            num_hap = num_haplotypes)
    
    sample_combo_likelihood <- NULL
    
    try(sample_combo_likelihood <- log(XNomial::xmulti(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$pProb), silent = TRUE)
    
    if (is.null(sample_combo_likelihood)) {
      sample_combo_likelihood <- log(XNomial::xmonte(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$observedProb)
    }
    
    combo_likelihood <- combo_likelihood + sample_combo_likelihood
    
  }
  
  ntrials_set <- 10000
  
  for (non_gene_sample in non_gene_samples) {
    
    exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = as.numeric(possible_strain_haplotype_combos),
                                                            strain_abun = as.numeric(phylotype_strains[non_gene_sample, ]),
                                                            num_hap = num_haplotypes)
    
    sample_combo_prob <- XNomial::xmonte(obs = sample_obs_abun[[non_gene_sample]], expr = exp_prob, ntrials = ntrials_set, detail = 0)$observedProb
    
    if (sample_combo_prob == 0) {
      sample_combo_prob <- 1 / ntrials_set
    }
    
    combo_likelihood <- combo_likelihood + log(sample_combo_prob)
    
  }
  
  return(-1 * combo_likelihood)
  
}




###res <- maxLik(logLik = combo_loglik,
###              start = c(10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10))





combo_loglik_gene_strain_intersect <- function(possible_strain_haplotype_combos_rep) {
  
  combo_likelihood <- 0
  
  for (gene_sample in gene_and_strain_samples) {
    
    exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = possible_strain_haplotype_combos_rep,
                                                            strain_abun = as.numeric(phylotype_strains[gene_sample, ]),
                                                            num_hap = num_haplotypes)
    
    sample_combo_likelihood <- NULL
    
    try(sample_combo_likelihood <- log(XNomial::xmulti(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$pProb), silent = TRUE)
    
    if (is.null(sample_combo_likelihood)) {
      sample_combo_likelihood <- log(XNomial::xmonte(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$observedProb)
    }
    
    combo_likelihood <- combo_likelihood + sample_combo_likelihood
    
  }
  
  ntrials_set <- 10000
  
  for (non_gene_sample in non_gene_samples) {
    
    exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = as.numeric(possible_strain_haplotype_combos_rep),
                                                            strain_abun = as.numeric(phylotype_strains[non_gene_sample, ]),
                                                            num_hap = num_haplotypes)
    
    sample_combo_prob <- XNomial::xmonte(obs = sample_obs_abun[[non_gene_sample]], expr = exp_prob, ntrials = ntrials_set, detail = 0)$observedProb
    
    if (sample_combo_prob == 0) {
      sample_combo_prob <- 1 / ntrials_set
    }
    
    combo_likelihood <- combo_likelihood + log(sample_combo_prob)
    
  }
  
  return(combo_likelihood)
  
}



combo_loglik_gene_strain_intersect_binom <- function(possible_strain_haplotype_combos_rep) {
  
  combo_likelihood <- 0
  
  for (gene_sample in gene_and_strain_samples) {
    
    exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = possible_strain_haplotype_combos_rep,
                                                            strain_abun = as.numeric(phylotype_strains[gene_sample, ]),
                                                            num_hap = num_haplotypes)
    
    sample_combo_likelihood <- NULL
    
    try(sample_combo_likelihood <- log(XNomial::xmulti(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$pProb), silent = TRUE)
    
    if (is.null(sample_combo_likelihood)) {
      sample_combo_likelihood <- log(XNomial::xmonte(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$observedProb)
    }
    
    combo_likelihood <- combo_likelihood + sample_combo_likelihood
    
  }
  
  for (non_gene_sample in non_gene_samples) {
    
    strains_w_gene <- which(as.numeric(possible_strain_haplotype_combos_rep) != num_haplotypes + 1)
    
    if (length(strains_w_gene) == 0) {
      next
    }
    
    exp_prob <- sum(phylotype_strains[non_gene_sample, strains_w_gene]) / sum(phylotype_strains[non_gene_sample, ])
    
    binom_test_p <- binom.test(x = 0, n = other_mean_variant_depth[[non_gene_sample]], p = exp_prob)$p.value
    
    if (binom_test_p < 1e-30) {
      binom_test_p <- 1e-30
    }
    
    combo_likelihood <- combo_likelihood + log(binom_test_p)
  }
  
  return(combo_likelihood)
  
}

combo_loglik_gene_strain_intersect_only <- function(possible_strain_haplotype_combos_rep) {
  
  combo_likelihood <- 0
  
  for (gene_sample in gene_and_strain_samples) {
    
    exp_prob <- compute_haplotype_exp_prob_based_on_strains(haplotype_dist = possible_strain_haplotype_combos_rep,
                                                            strain_abun = as.numeric(phylotype_strains[gene_sample, ]),
                                                            num_hap = num_haplotypes)
    
    sample_combo_likelihood <- NULL
    
    try(sample_combo_likelihood <- log(XNomial::xmulti(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$pProb), silent = TRUE)
    
    if (is.null(sample_combo_likelihood)) {
      sample_combo_likelihood <- log(XNomial::xmonte(obs = sample_obs_abun[[gene_sample]], expr = exp_prob, detail = 0)$observedProb)
    }
    
    combo_likelihood <- combo_likelihood + sample_combo_likelihood
    
  }
  
  return(combo_likelihood)
  
}


# Calculate aitchison distance between observed and expected haplotype abundances.
clr_transform_by_col <- function(in_df) {
  # Perform CLR transformation on table assuming samples are columns.
  return(data.frame(apply(in_df, 2, function(x){log(x) - mean(log(x))}), check.names=FALSE))
}

aitchison_dist_haplotypes <- function(strain_inferred_haplotypes) {
  
  haplotype_abun <- otu_tab[gene_and_strain_samples, ]
  haplotype_abun <- data.frame(sweep(haplotype_abun, 1, rowSums(haplotype_abun), '/'))
  colnames(haplotype_abun) <- gsub("^V", "H", colnames(haplotype_abun))
  
  phylotype_strains_raw <- phylotype_strains[gene_and_strain_samples, ]
  inferred_strain_link_haplotype_abun <- otu_tab[gene_and_strain_samples, ]
  inferred_strain_link_haplotype_abun[inferred_strain_link_haplotype_abun >= 0] <- 0.01
  for (strain_i in 1:length(strain_inferred_haplotypes)) {
    strain_haplotype <- strain_inferred_haplotypes[strain_i]
    
    if (strain_haplotype == num_haplotypes + 1) {
      next
    } else {
      inferred_strain_link_haplotype_abun[rownames(phylotype_strains_raw), strain_haplotype] <- phylotype_strains_raw[, strain_i]
    }
    
  }
  
  inferred_strain_link_haplotype_abun <- data.frame(sweep(inferred_strain_link_haplotype_abun, 1, rowSums(inferred_strain_link_haplotype_abun), '/'))
  colnames(inferred_strain_link_haplotype_abun) <- gsub("^V", "H", colnames(inferred_strain_link_haplotype_abun))
  
  haplotype_abun_clr <- clr_transform_by_col(t(haplotype_abun))
  inferred_strain_link_haplotype_abun_clr <-  clr_transform_by_col(t(inferred_strain_link_haplotype_abun))
  
  aitch_dist <- c()
  
  for (s in colnames(haplotype_abun_clr)) {
    aitch_dist <- c(aitch_dist, sqrt(sum((haplotype_abun_clr[, s] - inferred_strain_link_haplotype_abun_clr[, s])**2)))
  }
  
  return(mean(aitch_dist))
}

possible_strain_haplotype_combos_split <- as.list(data.frame(t(possible_strain_haplotype_combos)))


start_time <- Sys.time()
aitchison_dist_haplotypes_out <- mclapply(X = possible_strain_haplotype_combos_split, FUN = aitchison_dist_haplotypes, mc.cores = 10)
Sys.time() - start_time
aitchison_dist_haplotypes_out_vec <- unlist(aitchison_dist_haplotypes_out)


start_time <- Sys.time()
combo_loglik <- mclapply(X = possible_strain_haplotype_combos_split, FUN = combo_loglik_gene_strain_intersect_only, mc.cores = 10)
Sys.time() - start_time

combo_loglik_vec <- unlist(combo_loglik)



plot(aitchison_dist_haplotypes_out_vec, combo_loglik_vec)

combo_loglik_vec_AIC <- (-2 * combo_loglik_vec)

distinct_30_col <- c("#620045",
                     "#68ce5c",
                     "#72007a",
                     "#ace56d",
                     "#715bd1",
                     "#d3dc5f",
                     "#9d0181",
                     "#01b164",
                     "#ff7adc",
                     "#006a1e",
                     "#c98dff",
                     "#9e9b00",
                     "#014997",
                     "#ffb358",
                     "#e1a8ff",
                     "#3e7200",
                     "#f54c98",
                     "#006a3a",
                     "#c91f3e",
                     "#c6de84",
                     "#a7004f",
                     "#ffcb7e",
                     "#9668a7",
                     "#c25a04",
                     "#ffaade",
                     "#a05b00",
                     "#ba6794",
                     "#662c00",
                     "#ff9580",
                     "#b60045")

haplotype_cols <- distinct_30_col[1:num_haplotypes]

NA_col <- "#abaaaa"


strain_inferred_haplotypes <- possible_strain_haplotype_combos_split[[388]]
#strain_inferred_haplotypes <- possible_strain_haplotype_combos_split[[which.min(aitchison_dist_haplotypes_out)]]

# Make plot showing the haplotype relative abundances straight from StrainFinder
# (including for samples where strains are present, but haplotypes aren't).

haplotype_abun <- otu_tab[gene_and_strain_samples, ] - 0.01
haplotype_abun <- data.frame(sweep(haplotype_abun, 1, rowSums(haplotype_abun), '/'))
haplotype_abun[is.na(haplotype_abun)] <- 0
colnames(haplotype_abun) <- gsub("^V", "H", colnames(haplotype_abun))
haplotype_abun$Sample <- factor(rownames(haplotype_abun), levels = rownames(haplotype_abun))
haplotype_abun <- melt(haplotype_abun,
                       id.vars = "Sample",
                       variable.name = "Haplotype",
                       value.name = "rel_abun")

haplotype_abun_plot <- ggplot(data = haplotype_abun, aes(y = Sample, x = rel_abun, fill = Haplotype)) +
  geom_bar(position="stack", stat="identity") +
  #scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(values = haplotype_cols) +
  ylab("Sample") +
  xlab("Relative abundance") +
  theme(legend.position = "none")



# Then make plot after computing the haplotype abundances based on the abundances of strains 
phylotype_strains_raw <- phylotype_strains - 0.01
phylotype_strains_raw <- phylotype_strains_raw[gene_and_strain_samples, ]

inferred_strain_link_haplotype_abun <- otu_tab[gene_and_strain_samples, ]
inferred_strain_link_haplotype_abun[inferred_strain_link_haplotype_abun > 0] <- 0

for (strain_i in 1:length(strain_inferred_haplotypes)) {
  strain_haplotype <- strain_inferred_haplotypes[strain_i]
  
  if (strain_haplotype == num_haplotypes + 1) {
    next
  } else {
    inferred_strain_link_haplotype_abun[rownames(phylotype_strains_raw), strain_haplotype] <- phylotype_strains_raw[, strain_i]
  }
  
}

inferred_strain_link_haplotype_abun <- data.frame(sweep(inferred_strain_link_haplotype_abun, 1, rowSums(inferred_strain_link_haplotype_abun), '/'))
inferred_strain_link_haplotype_abun[is.na(inferred_strain_link_haplotype_abun)] <- 0
colnames(inferred_strain_link_haplotype_abun) <- gsub("^V", "H", colnames(inferred_strain_link_haplotype_abun))
inferred_strain_link_haplotype_abun$Sample <- factor(rownames(inferred_strain_link_haplotype_abun), levels = rownames(inferred_strain_link_haplotype_abun))
inferred_strain_link_haplotype_abun <- melt(inferred_strain_link_haplotype_abun,
                                            id.vars = "Sample",
                                            variable.name = "Haplotype",
                                            value.name = "rel_abun")

inferred_strain_link_haplotype_abun_plot <- ggplot(data = inferred_strain_link_haplotype_abun, aes(y = Sample, x = rel_abun, fill = Haplotype)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = haplotype_cols) +
  ylab("") +
  xlab("Relative abundance")


legend <- get_legend(inferred_strain_link_haplotype_abun_plot)

# and replot suppressing the legend
inferred_strain_link_haplotype_abun_plot <- inferred_strain_link_haplotype_abun_plot + theme(legend.position='none')


# Make basic plot showing how haplotypes were inferred to be distributed across strains.

strain_haplotype_links <- data.frame(Strain = 1:length(strain_inferred_haplotypes),
                                     Haplotype = strain_inferred_haplotypes)

if (length(which(strain_haplotype_links$Haplotype == num_haplotypes + 1)) > 0) {
  strain_haplotype_links[which(strain_haplotype_links$Haplotype == num_haplotypes + 1), "Haplotype"] <- NA
}

strain_haplotype_links$Strain <- factor(strain_haplotype_links$Strain, levels = rev(strain_haplotype_links$Strain))
strain_haplotype_links$Haplotype <- colnames(otu_tab)[strain_haplotype_links$Haplotype]
present_haplotypes_col <- haplotype_cols[which(colnames(otu_tab) %in% strain_haplotype_links$Haplotype)]
strain_haplotype_links$Haplotype <- gsub("V", "H", strain_haplotype_links$Haplotype)


strain_haplotype_links_plot <- ggplot(data = strain_haplotype_links, aes(y = Strain, x = 1, fill = Haplotype)) +
  geom_tile() +
  scale_fill_manual(values = present_haplotypes_col) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# now make combined plot:
ggdraw(plot_grid(plot_grid(strain_haplotype_links_plot, haplotype_abun_plot, inferred_strain_link_haplotype_abun_plot, ncol=3, rel_widths = c(0.2, 0.4, 0.4)),
                 plot_grid(NULL, legend, ncol=1), align = "v",
                 rel_widths=c(1, 0.2)))

