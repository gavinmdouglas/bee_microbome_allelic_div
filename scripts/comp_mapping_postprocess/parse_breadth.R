### Parse breadth of coverage from bedGraph files and add in useful summary info.

rm(list = ls(all.names = TRUE))

library(plyr)

read_in_breadth_files <- function(in_path, pattern, roary_formatted_pangenome) {
  
  in_path <- gsub("/$", "", in_path)
  
  # Read in all breadth files.
  input_breadth_files <- list.files(path = in_path, full.names = TRUE, pattern = pattern)
  
  input_breadth_by_sample <- lapply(input_breadth_files, function(x) { read.table(x, header = FALSE, sep = "\t", row.names = 1, stringsAsFactors = FALSE) })
  input_samples <- gsub(paste(in_path, "/", sep = ""), "", input_breadth_files)
  input_samples <- gsub(pattern, "", input_samples)
  names(input_breadth_by_sample) <- input_samples
  
  input_breadth_by_sample[[1]]$sample <- input_samples[1]
  input_breadth_by_sample[[1]]$gene <- rownames(input_breadth_by_sample[[1]])
  rownames(input_breadth_by_sample[[1]]) <- NULL
  input_breadth_combined <- input_breadth_by_sample[[1]]
  
  for (i in 2:length(input_samples)) {
    input_breadth_by_sample[[i]]$sample <- input_samples[i]
    input_breadth_by_sample[[i]]$gene <- rownames(input_breadth_by_sample[[i]])
    rownames(input_breadth_by_sample[[i]]) <- NULL
    input_breadth_combined <- rbind(input_breadth_combined,  input_breadth_by_sample[[i]])
  }
  
  input_breadth_combined <- input_breadth_combined[, c("sample", "gene", "V7")]
  colnames(input_breadth_combined)[3] <- "breadth"
  
  input_unique_genes <- unique(input_breadth_combined$gene)
  
  # Get per-gene summary across samples.
  input_breadth_summary <- data.frame(matrix(NA,
                                             nrow = length(input_unique_genes),
                                             ncol = 13))
  colnames(input_breadth_summary) <- c("ref_num", "mean", "median", "max", "min", "sd", "cv", "num_at_least_0.1",
                                       "num_at_least_0.25", "num_at_least_0.5", "num_at_least_0.9", "num_1", "num_0")
  rownames(input_breadth_summary) <- input_unique_genes
  
  for (input_gene in input_unique_genes) {
    input_breadth_summary[input_gene, "ref_num"] <- roary_formatted_pangenome[input_gene, "No..isolates"]
    input_gene_breadth_values <- input_breadth_combined[which(input_breadth_combined$gene == input_gene), "breadth"]
    
    input_breadth_summary[input_gene, "mean"] <- mean(input_gene_breadth_values)
    input_breadth_summary[input_gene, "median"] <- median(input_gene_breadth_values)
    input_breadth_summary[input_gene, "max"] <- max(input_gene_breadth_values)
    input_breadth_summary[input_gene, "min"] <- min(input_gene_breadth_values)
    input_breadth_summary[input_gene, "sd"] <- sd(input_gene_breadth_values)
    input_breadth_summary[input_gene, "cv"] <- input_breadth_summary[input_gene, "sd"] / input_breadth_summary[input_gene, "mean"]
    input_breadth_summary[input_gene, "num_at_least_0.1"] <- length(which(input_gene_breadth_values >= 0.1))
    input_breadth_summary[input_gene, "num_at_least_0.25"] <- length(which(input_gene_breadth_values >= 0.25))
    input_breadth_summary[input_gene, "num_at_least_0.5"] <- length(which(input_gene_breadth_values >= 0.5))
    input_breadth_summary[input_gene, "num_at_least_0.9"] <- length(which(input_gene_breadth_values >= 0.9))
    input_breadth_summary[input_gene, "num_1"] <- length(which(input_gene_breadth_values == 1))
    input_breadth_summary[input_gene, "num_0"] <- length(which(input_gene_breadth_values == 0))
  }
  
  # Create rounded version of metrics.
  input_breadth_summary$mean_rounded <- round(input_breadth_summary$mean, digits = 1)
  input_breadth_summary$median_rounded <- round(input_breadth_summary$median, digits = 1)
  input_breadth_summary$max_rounded <- round(input_breadth_summary$max, digits = 1)
  input_breadth_summary$min_rounded <- round(input_breadth_summary$min, digits = 1)
  input_breadth_summary$sd_rounded <- round(input_breadth_summary$sd, digits = 1)
  input_breadth_summary$cv_rounded <- round(input_breadth_summary$cv, digits = 1)
  
  return(list(breadth_by_sample = input_breadth_by_sample,
              breadth_all_samples = input_breadth_combined,
              breadth_summary = input_breadth_summary,
              unique_genes = input_unique_genes))
  
}

combined_panaroo <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/highqual_genomes_panaroo/combined_panaroo_annot.rds")

combined_panaroo[which(is.na(combined_panaroo$No..isolates)), "No..isolates"] <- 1

breadth_output <- read_in_breadth_files(in_path = "/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/comp_mapping/coverage_breadth",
                                        pattern = ".cov.bedGraph.gz",
                                        roary_formatted_pangenome = combined_panaroo)


breadth_output$breadth_summary$species <- NA
all_species <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/species.txt", stringsAsFactors = FALSE)$V1
for (sp in all_species) {
  breadth_output$breadth_summary[grep(sp, rownames(breadth_output$breadth_summary)), "species"] <- sp
}


panaroo_only_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/panaroo_only_potential_core.rds")
checkm_panaroo_potential_core <- readRDS("/data1/gdouglas/projects/honey_bee/ref_genome_pangenomes/species/core_genes/checkm_panaroo_passed_potential_core.rds")
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
  breadth_by_sample_clean[[sample]] <- breadth_output$breadth_by_sample[[sample]][, c("gene", "V7"), drop = FALSE]
  colnames(breadth_by_sample_clean[[sample]]) <- c("gene", sample)
}


# Call genes as present based on breadth of coverage of at least 0.5.
all_present <- join_all(breadth_by_sample_clean, by="gene", type='left')
rownames(all_present) <- all_present$gene
all_present <- all_present[, -which(colnames(all_present) == "gene")]
all_present[all_present < 0.5] <- 0
all_present[all_present >= 0.5] <- 1



saveRDS(object = all_present,
        file = "/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/comp_mapping/summary/gene_presence_0.5_breadth.rds")

saveRDS(object = breadth_output,
        file = "/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/comp_mapping/summary/gene_breadth_breakdown.rds")

