rm(list = ls(all.names = TRUE))

# Create mapfile of Prokka genes to each genome.

all_species <- read.table('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/final_species_names.txt.gz',
                          header = FALSE, stringsAsFactors = FALSE)$V1

# Treat Apilactobacillus_apinorum separately.
Aa_panaroo_file <- '/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_panaroo/Apilactobacillus_apinorum/gene_presence_absence.csv.gz'
Aa_gene_to_genome <- read.table(file = Aa_panaroo_file, header = TRUE, sep = ',', stringsAsFactors = FALSE, quote = '', comment.char = '')
Apilactobacillus_apinorum_out <- data.frame(gene=Aa_gene_to_genome$GCA_001281175.1,
                                            genome='GCA_001281175.1')

write.table(x = Apilactobacillus_apinorum_out,
            file = '/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_to_genome_by_species/Apilactobacillus_apinorum.tsv',
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')

for (sp in all_species) {
  
  if (sp == 'Apilactobacillus_apinorum') { next }

  sp_gene_info_file <- paste('/data1/gdouglas/projects/honey_bee/ref_genomes/highqual_genomes_panaroo/', sp, '/gene_data.csv.gz', sep = '')
  sp_gene_to_genome <- read.table(file = sp_gene_info_file, header = TRUE, sep = ',', stringsAsFactors = FALSE, quote = '', comment.char = '')

  output <- data.frame(gene=sp_gene_to_genome$annotation_id,
                       genome=sp_gene_to_genome$gff_file)
  
  write.table(x = output,
              file = paste('/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/gene_to_genome_by_species/', sp, '.tsv', sep = ''),
              col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t')
}
