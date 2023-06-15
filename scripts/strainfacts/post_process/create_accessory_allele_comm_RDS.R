rm(list = ls(all.names = TRUE))

# Create gene allele rel. abun. tables converted to lists of RDS objects that are quicker to read in and use.
library(parallel)

setwd('/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/comp_mapping/strainfacts_running/accessory/output/')

preprocess_comm <- function(gene) {

  gene_samples <- read.table(paste('/data1/gdouglas/projects/honey_bee/large_files_to_backup/strainfacts/prepped_input/prepped_accessory_input/samples/', gene, '_samples.tsv', sep = ''),
                             sep = "\t", row.names = 1, header = FALSE)

  allele_freq <- read.table(paste('comm/', gene, ".comm.tsv.gz", sep = ""),
                               header = TRUE, sep = "\t")
  
  if (length(unique(allele_freq$sample)) != nrow(gene_samples)) {
    stop('Mismatch in expected number of samples and those in sample mapfile!') 
  }
  
  allele_freq$sample <- as.character(allele_freq$sample)
  
  if (length(which(! allele_freq$sample %in% rownames(gene_samples))) > 0) {
    stop('Error - sample missing in mapfile.') 
  }
  
  allele_freq$sample <- gene_samples[allele_freq$sample, 1]
  
  orig_samples <- sort(unique(allele_freq$sample))
  
  # Only keep strain calls that were at least >= 1%
  allele_freq <- allele_freq[which(allele_freq$community >= 0.01), ]
  
  if (! identical(orig_samples, sort(unique(allele_freq$sample)))) {
     stop("ERROR - samples dropped after excluding low abun strains.")
  }
  
  # Total-sum scale by sample remaining rel. abun.
  sample_sums <- aggregate(x = community ~ sample,
                           data = allele_freq,
                           FUN = sum)
  rownames(sample_sums) <- sample_sums$sample
  allele_freq$community <- (allele_freq$community / sample_sums[allele_freq$sample, "community"]) * 100
  
  return(allele_freq)

}

all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          header = FALSE, stringsAsFactors = FALSE)$V1

allele_relabun <- list()

all_fit <- list.files("fit/", pattern = ".fit$")
all_genes <- gsub(".fit$", "", all_fit)

for (sp in all_species) {
  
  print(sp)
  
  sp_genes <- grep(sp, all_genes, value = TRUE)
  
  if (length(sp_genes) == 0) { next }
  
  allele_relabun[[sp]] <- parallel::mclapply(X = sp_genes, FUN = preprocess_comm, mc.cores = 40)

  names(allele_relabun[[sp]]) <- sp_genes

}

saveRDS(object = allele_relabun,
        file = "/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/RDS/strainfacts_accessory_allele_relabun.rds")
