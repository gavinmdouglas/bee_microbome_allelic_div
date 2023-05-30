# Parse breadth of coverage from bedGraph files and add in useful summary info.
# Define genes as present based on calls of >= 50% coverage.

rm(list = ls(all.names = TRUE))

library(parallel)
library(plyr)

read_in_breadth_files <- function(in_path, pattern, roary_formatted_pangenome, num_cores = 1) {

  in_path <- gsub("/$", "", in_path)
  
  # Read in all breadth files.
  input_breadth_files <- list.files(path = in_path, full.names = TRUE, pattern = pattern)
  
  input_breadth_by_sample <- parallel::mclapply(input_breadth_files,
                                                function(x) { read.table(x, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE) },
                                                mc.cores = num_cores)
  
  input_samples <- gsub(paste(in_path, "/", sep = ""), "", input_breadth_files)
  input_samples <- gsub(pattern, "", input_samples)
  names(input_breadth_by_sample) <- input_samples
  
  input_breadth_by_sample_added <- parallel::mclapply(input_samples,
                     function(x) {
                                   working_df <-  input_breadth_by_sample[[x]]
                                   working_df$sample <- x
                                   working_df$gene <- rownames(working_df)
                                   rownames(working_df) <- NULL
                                   working_df <- working_df[, c("sample", "gene", "V7")]
                                   colnames(working_df)[3] <- "breadth"
                                   return(working_df)
                                },
                    mc.cores = num_cores)
  
  input_breadth_by_sample <- input_breadth_by_sample_added
  names(input_breadth_by_sample) <- input_samples
  rm(input_breadth_by_sample_added)
  
  input_breadth_combined <- data.table::rbindlist(input_breadth_by_sample)
  
  input_breadth_combined <- data.frame(input_breadth_combined)
  
  input_unique_genes <- unique(input_breadth_combined$gene)
  
  input_unique_gene_metrics <- parallel::mclapply(input_unique_genes,
                                                  function(gene) {
                                                    
                                                    gene_breadth_values <- as.numeric(input_breadth_combined[which(input_breadth_combined$gene == gene), "breadth"])
                                                    
                                                    return(
                                                      
                                                      data.frame(matrix(c(roary_formatted_pangenome[gene, "No..isolates"],
                                                                          mean(gene_breadth_values),
                                                                          median(gene_breadth_values),
                                                                          max(gene_breadth_values),
                                                                          min(gene_breadth_values),
                                                                          sd(gene_breadth_values),
                                                                          sd(gene_breadth_values) / mean(gene_breadth_values),
                                                                          length(which(gene_breadth_values >= 0.1)),
                                                                          length(which(gene_breadth_values >= 0.25)),
                                                                          length(which(gene_breadth_values >= 0.5)),
                                                                          length(which(gene_breadth_values >= 0.9)),
                                                                          length(which(gene_breadth_values == 1)),
                                                                          length(which(gene_breadth_values == 0))),
                                                                         nrow = 1, ncol = 13)
                                                    ))
                                                  },
                                                  mc.cores = num_cores)
  
  # Get per-gene summary across samples.
  input_breadth_summary <- data.table::rbindlist(input_unique_gene_metrics)
  input_breadth_summary <- data.frame(input_breadth_summary)
  colnames(input_breadth_summary) <- c("ref_num", "mean", "median", "max", "min", "sd", "cv", "num_at_least_0.1",
                                       "num_at_least_0.25", "num_at_least_0.5", "num_at_least_0.9", "num_1", "num_0")
  rownames(input_breadth_summary) <- input_unique_genes
  
  return(list(breadth_by_sample = input_breadth_by_sample,
              breadth_all_samples = input_breadth_combined,
              breadth_summary = input_breadth_summary,
              unique_genes = input_unique_genes))
  
}

combined_panaroo <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/combined_panaroo_annot.rds")
combined_panaroo[which(is.na(combined_panaroo$No..isolates)), "No..isolates"] <- 1

breadth_output <- read_in_breadth_files(in_path = "/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/coverage_breadth",
                                        pattern = ".cov.bedGraph.gz",
                                        roary_formatted_pangenome = combined_panaroo,
                                        num_cores = 50)

breadth_output$breadth_summary$species <- NA
all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          stringsAsFactors = FALSE)$V1
for (sp in all_species) {
  breadth_output$breadth_summary[grep(sp, rownames(breadth_output$breadth_summary)), "species"] <- sp
}

panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/core_genes/RDS/checkm_panaroo_passed_potential_core.rds")
core_genes <- c(panaroo_only_potential_core, checkm_panaroo_potential_core)

breadth_output$breadth_summary$gene_type <- "Non-core"
for (sp in all_species) {
  breadth_output$breadth_summary[core_genes[[sp]], "gene_type"] <- "Core"
}

# Remove core genes that shouldn't have been added in
# (because they were removed after trimming genes in bed due to being too short)
breadth_output$breadth_summary <- breadth_output$breadth_summary[which(rowSums(is.na(breadth_output$breadth_summary)) != 19), ]

breadth_by_sample_clean <- list()

for (sample in names(breadth_output$breadth_by_sample)) {
  breadth_by_sample_clean[[sample]] <- breadth_output$breadth_by_sample[[sample]][, c('gene', 'breadth'), drop = FALSE]
  colnames(breadth_by_sample_clean[[sample]]) <- c('gene', sample)
}


# Call genes as present based on breadth of coverage of at least 0.5.
all_present <- plyr::join_all(breadth_by_sample_clean, by="gene", type='left')
rownames(all_present) <- all_present$gene
all_present <- all_present[, -which(colnames(all_present) == "gene")]

# Also keep track of the matrix of raw breadths, in case you want to explore different cut-offs more later.
all_breadth <- all_present

all_present[all_present < 0.5] <- 0
all_present[all_present >= 0.5] <- 1

write.table(x = all_present,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_presence_0.5_breadth.tsv',
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = NA)

write.table(x = all_breadth,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/comp_mapping/summary/gene_breadth.tsv',
            quote = FALSE, sep = '\t', row.names = TRUE, col.names = NA)

saveRDS(object = breadth_output,
        file = '/data1/gdouglas/projects/honey_bee/large_files_to_backup/mgs_datasets/comp_mapping/summary/gene_breadth_breakdown.rds')
