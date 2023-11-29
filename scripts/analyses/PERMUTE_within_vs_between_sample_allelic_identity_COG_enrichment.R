rm(list = ls(all.names = TRUE))

source('/home/gdouglas/scripts/bee_microbome_allelic_div/scripts/functions.R')

results <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/mgs_datasets/strainfacts/accessory/within_sample_allele_permutation_results.tsv.gz',
                      header = TRUE, sep = '\t', stringsAsFactors = FALSE)

results$species <- cleanup_species_names(results$species)

results <- results[which(results$num_alleles >= 5), ]
results <- results[which(results$num_samples >= 5), ]
results <- results[which(results$mean_alleles_per_sample >= 1.5), ]

results$p_lower_BH <- p.adjust(results$p_lower, 'BH')
results$p_higher_BH <- p.adjust(results$p_higher, 'BH')

COG_category_info <- read.table(file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/combined_pangenome/eggnog_mapper_annot.tsv.gz',
                                header = TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1, quote = "", comment.char = "")
COG_category_info[which(COG_category_info$COG_category == ''), 'COG_category'] <- '-'

COG_gene_to_category <- read.table("/data1/gdouglas/db/COG_definitions/cog-20.to_category.tsv", header = FALSE, sep = "\t")
COG_category_to_COG <- list()
for (category in unique(COG_gene_to_category$V2)) {
  COG_category_to_COG[[category]] <- COG_gene_to_category[which(COG_gene_to_category$V2 == category), "V1"]
}

COG_descrip <- read.table('/data1/gdouglas/db/COG_definitions/cog-20.def.tab',
                          header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)
COG_category_descrip <- read.table('/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv',
                                   header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

COG_category_descrip$clean <- paste(rownames(COG_category_descrip), COG_category_descrip$V2, sep = ' - ')
COG_category_descrip['-', 'clean'] <- 'No COG annotation'

datasets <- sort(unique(results$dataset))
all_species <- sort(unique(results$species))
categories_to_ignore <- c('A', 'B', 'Y', 'Z')

clean_COG_vec <- function(in_vec) {
  in_vec <- in_vec[which(! is.na(in_vec))]
  in_vec <- in_vec[which(in_vec != '')]
  
  multi_cog_elements <- grep(',', in_vec)
  if (length(multi_cog_elements) > 0) {
    multi_subset <- in_vec[multi_cog_elements]
    in_vec <- in_vec[-multi_cog_elements]
    
    for (multi_cog_element in multi_subset) {
      in_vec <- c(in_vec, strsplit(multi_cog_element, ',')[[1]])
    }
  }

  return(in_vec)
}

raw_out <- list()

for (sp in all_species) {
  
  for (d in datasets) {
    
    results_subset <- results[which(results$species == sp & results$dataset == d), ]
    
    sig_genes <- results_subset[which(results_subset$p_higher_BH < 0.25), 'gene']
    nonsig_genes <- results_subset[which(results_subset$p_higher_BH >= 0.25), 'gene']

  
    if (length(sig_genes) < 10 || length(nonsig_genes) < 10) { next }
  
    sig_raw_COG <- COG_category_info[sig_genes, "all_COG"]
    nonsig_raw_COG <- COG_category_info[nonsig_genes, "all_COG"]
      
    sig_COGs <- clean_COG_vec(sig_raw_COG)
    nonsig_COGs <- clean_COG_vec(nonsig_raw_COG)

    enrichment_output <- identify_enriched_categories(genes = sig_COGs,
                                                   background = nonsig_COGs,
                                                   gene_to_category_map = COG_category_to_COG,
                                                   min_category_count = 10,
                                                   to_ignore = categories_to_ignore)
    
    if (nrow(enrichment_output) == 0) { next }
    
    enrichment_output$species <- sp
    enrichment_output$datset <- d
    
    raw_out[[paste(sp, d)]] <- enrichment_output
    
 }

}

tmp <- do.call(rbind, raw_out)
tmp$descrip <- COG_category_descrip[tmp$category, 'clean']

tmp_sig <- tmp[which(tmp$fdr < 0.25), ]

write.table(x = tmp_sig, file = '/home/gdouglas/tmp/higher_indentity_sig_COG_category.tsv', col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
