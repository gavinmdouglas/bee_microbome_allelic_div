rm(list = ls(all.names = TRUE))

# ANI values to report in text.
fastani_out <- read.table("/data1/gdouglas/projects/bee_microbiome_zenodo/ref_genomes/CheckM-filt_FastANI_output.tsv.gz",
                          header = FALSE, sep = "\t", stringsAsFactors = FALSE)

fastani_out$V1 <- gsub("\\.fna$", "", fastani_out$V1)
fastani_out$V2 <- gsub("\\.fna$", "", fastani_out$V2)
fastani_out$V1 <- gsub("../../highqual_genomes_draft/", "", fastani_out$V1)
fastani_out$V2 <- gsub("../../highqual_genomes_draft/", "", fastani_out$V2)

colnames(fastani_out) <- c("file1", "file2", "ani", "num_aligned", "num_fragments")

fastani_out$sp1 <- gsub("/.*$", "", fastani_out$file1)
fastani_out$id1 <- gsub("^.*/", "", fastani_out$file1)

fastani_out$sp2 <- gsub("/.*$", "", fastani_out$file2)
fastani_out$id2 <- gsub("^.*/", "", fastani_out$file2)

# Remove instances comparisons between the same genome.
fastani_out <- fastani_out[which(fastani_out$id1 != fastani_out$id2), ]

# Also, only keep first instance of every pairwise comparison.
fastani_out$id_combined <- NA
for (i in 1:nrow(fastani_out)) {
  id_combined <- sort(c(fastani_out[i, "id1"], fastani_out[i, "id2"]))
  fastani_out$id_combined[i] <- paste(id_combined, collapse = " ")
}
fastani_out <- fastani_out[which(! duplicated(fastani_out$id_combined)), ]

# Function to get table subset to species comparisons *only* within and between the specified set.
subset_fastani <- function(fastani_df, species_subset) {
 
  rows_to_keep <- as.numeric()
  
  for (i in 1:nrow(fastani_df)) {
    if ((fastani_df[i, 'sp1'] %in% species_subset) && (fastani_df[i, 'sp2'] %in% species_subset)) {
      rows_to_keep <- c(rows_to_keep, i)
    }
  }
  
  return(fastani_df[rows_to_keep, , drop = FALSE])
}

# Six Lactobacillus genomes clustering with L. kullabergensis and L. kimbladii.
kullabergensis_kimbladii_and_six_others <- subset_fastani(fastani_out, c('Lactobacillus_kullabergensis', 'Lactobacillus_kimbladii', 'Lactobacillus_sp'))

## Only keep rows comparing between these two species and Lactobacillus sp. 
kullabergensis_kimbladii_vs_six_others_i <- c(which((kullabergensis_kimbladii_and_six_others$sp1 == 'Lactobacillus_sp') & (kullabergensis_kimbladii_and_six_others$sp2 != 'Lactobacillus_sp')),
                                              which((kullabergensis_kimbladii_and_six_others$sp1 != 'Lactobacillus_sp') & (kullabergensis_kimbladii_and_six_others$sp2 == 'Lactobacillus_sp')))
kullabergensis_kimbladii_vs_six_others <- kullabergensis_kimbladii_and_six_others[kullabergensis_kimbladii_vs_six_others_i, ]

## Get max pairwise ANI for each of these six vs either L. kullabergensis or L. kimbladii.
kullabergensis_kimbladii_vs_six_others$Lactobacillus_sp_id <- kullabergensis_kimbladii_vs_six_others$id1
kullabergensis_kimbladii_vs_six_others$Lactobacillus_sp_id[which(kullabergensis_kimbladii_vs_six_others$sp2 == 'Lactobacillus_sp')] <- kullabergensis_kimbladii_vs_six_others[which(kullabergensis_kimbladii_vs_six_others$sp2 == 'Lactobacillus_sp'), 'id2']

Lactobacillus_sp_to_exclude <- c("GCA_000760615.1", "GCF_016101625.1", "GCA_003693025.1",
                                 "GCF_016100885.1", "GCF_016100935.1", "GCF_016100975.1")
kullabergensis_kimbladii_vs_six_others <- kullabergensis_kimbladii_vs_six_others[which(kullabergensis_kimbladii_vs_six_others$Lactobacillus_sp_id %in% Lactobacillus_sp_to_exclude), ]

round(range(aggregate(ani ~ Lactobacillus_sp_id, data = kullabergensis_kimbladii_vs_six_others, FUN = max)$ani), 1)


# Bombella
Bombella_ani <- subset_fastani(fastani_out, c('Saccharibacter_sp', 'Bombella_sp', 'Bombella_apis', 'Parasaccharibacter_apium'))

## Separate into two clusters.
Bombella_sp_genomes <- c('GCF_002592045.1', 'GCF_009725755.1', 'GCF_009725845.1')
Bombella_sp_genome_only_i <- which((Bombella_ani$id1 %in% Bombella_sp_genomes) & (Bombella_ani$id2 %in% Bombella_sp_genomes))
Bombella_sp_genome_excluded_i <- which((! Bombella_ani$id1 %in% Bombella_sp_genomes) & (! Bombella_ani$id2 %in% Bombella_sp_genomes))

Bombella_sp_comparisons <- Bombella_ani[Bombella_sp_genome_only_i, ]
Bombella_apis_comparisons <- Bombella_ani[Bombella_sp_genome_excluded_i, ]

Bombella_apis_comparisons_vs_others_i <- c(which(Bombella_apis_comparisons$sp1 == 'Bombella_apis' & Bombella_apis_comparisons$sp2 != 'Bombella_apis'),
                                           which(Bombella_apis_comparisons$sp1 != 'Bombella_apis' & Bombella_apis_comparisons$sp2 == 'Bombella_apis'))
Bombella_apis_comparisons_vs_others <- Bombella_apis_comparisons[Bombella_apis_comparisons_vs_others_i, ]
Bombella_apis_comparisons_vs_others$other_id <- Bombella_apis_comparisons_vs_others$id1
Bombella_apis_comparisons_vs_others$other_id[which(Bombella_apis_comparisons_vs_others$sp1 == 'Bombella_apis')] <- Bombella_apis_comparisons_vs_others[which(Bombella_apis_comparisons_vs_others$sp1 == 'Bombella_apis'), 'id2']

round(range(aggregate(ani ~ other_id, data = Bombella_apis_comparisons_vs_others, FUN = max)$ani), 1)

round(range(Bombella_sp_comparisons$ani), 1)


# GCF_016101275.1 vs Apilactobacillus kunkeei.
kunkeei_vs_additional_i <- c(which(fastani_out$sp1 == 'Apilactobacillus_kunkeei' & fastani_out$id2 == 'GCF_016101275.1'),
                             which(fastani_out$sp2 == 'Apilactobacillus_kunkeei' & fastani_out$id1 == 'GCF_016101275.1'))

round(range(fastani_out[kunkeei_vs_additional_i, 'ani']), 1)


# B. coryneforme vs B. indicum 
B_coryneforme_vs_indicum_i <- c(which(fastani_out$sp1 == 'Bifidobacterium_coryneforme' & fastani_out$sp2 == 'Bifidobacterium_indicum'),
                                which(fastani_out$sp1 == 'Bifidobacterium_indicum' & fastani_out$sp2 == 'Bifidobacterium_coryneforme'))
round(fastani_out[B_coryneforme_vs_indicum_i, 'ani'], 1)

# Within-B. coryneforme ANI 
B_coryneforme_within_i <- which(fastani_out$sp1 == 'Bifidobacterium_coryneforme' & fastani_out$sp2 == 'Bifidobacterium_coryneforme')
round(fastani_out[B_coryneforme_within_i, 'ani'], 1)


# Bombilactobacillus vs Lactobacillus sp. hits
Bombilactobacillus_mellis_lactobacilli <- c('GCF_016100965.1', 'GCF_016100985.1', 'GCF_016101055.1', 'GCF_016101025.1', 'GCF_016101075.1', 'GCF_016101045.1')
Bombilactobacillus_mellis_vs_lactobacilli_i <- c(intersect(grep('Bombilactobacillus_mellis', fastani_out$sp1), which(fastani_out$id2 %in% Bombilactobacillus_mellis_lactobacilli)),
                                          intersect(grep('Bombilactobacillus_mellis', fastani_out$sp2), which(fastani_out$id1 %in% Bombilactobacillus_mellis_lactobacilli)))
Bombilactobacillus_mellis_vs_lactobacilli <- fastani_out[Bombilactobacillus_mellis_vs_lactobacilli_i, ]
Bombilactobacillus_mellis_vs_lactobacilli$Lactobacillus_sp_id <- Bombilactobacillus_mellis_vs_lactobacilli$id1
Bombilactobacillus_mellis_vs_lactobacilli$Lactobacillus_sp_id[which(Bombilactobacillus_mellis_vs_lactobacilli$sp2 == 'Lactobacillus_sp')] <- Bombilactobacillus_mellis_vs_lactobacilli[which(Bombilactobacillus_mellis_vs_lactobacilli$sp2 == 'Lactobacillus_sp'), 'id2']
round(range(aggregate(ani ~ Lactobacillus_sp_id, data = Bombilactobacillus_mellis_vs_lactobacilli, FUN = max)$ani), 1)


Bombilactobacillus_mellifer_lactobacilli <- c('GCF_016102125.1')
Bombilactobacillus_mellifer_vs_lactobacilli_i <- c(intersect(grep('Bombilactobacillus_mellifer', fastani_out$sp1), which(fastani_out$id2 %in% Bombilactobacillus_mellifer_lactobacilli)),
                                                 intersect(grep('Bombilactobacillus_mellifer', fastani_out$sp2), which(fastani_out$id1 %in% Bombilactobacillus_mellifer_lactobacilli)))
Bombilactobacillus_mellifer_vs_lactobacilli <- fastani_out[Bombilactobacillus_mellifer_vs_lactobacilli_i, ]
round(Bombilactobacillus_mellifer_vs_lactobacilli$ani, 1)


# GCA_001723875.1 and GCA_007559165.1 vs Gilliamella apis
apis_recodes <- c('GCA_001723875.1', 'GCA_007559165.1')
Gilliamella_apis_vs_recodes_i <- c(which((fastani_out$sp1 == 'Gilliamella_apis') & (fastani_out$id2 %in% apis_recodes)),
                                   which((fastani_out$id1 %in% apis_recodes) & (fastani_out$sp2 == 'Gilliamella_apis')))
aggregate(ani ~ id1, data = fastani_out[Gilliamella_apis_vs_recodes_i, ], FUN = max)


# GCA_003202915.1 and GCA_003202655.1 vs Gilliamella apicola
apicola_recodes <- c('GCA_003202915.1', 'GCA_003202655.1')
Gilliamella_apicola_vs_recodes_i <- unique(c(which((fastani_out$sp1 == 'Gilliamella_apicola') & (fastani_out$id2 %in% apicola_recodes)),
                                             which((fastani_out$id1 %in% apicola_recodes) & (fastani_out$sp2 == 'Gilliamella_apicola'))))

Gilliamella_apicola_vs_recodes <- fastani_out[Gilliamella_apicola_vs_recodes_i, ]

within_apicola_recodes_i <- which(Gilliamella_apicola_vs_recodes$id1 %in% apicola_recodes & Gilliamella_apicola_vs_recodes$id2 %in% apicola_recodes)
within_recodes <- Gilliamella_apicola_vs_recodes[within_apicola_recodes_i, , drop = FALSE]
round(within_recodes$ani, 1)

Gilliamella_apicola_vs_recodes_only <- Gilliamella_apicola_vs_recodes[-within_apicola_recodes_i, ]
Gilliamella_apicola_vs_recodes_only$outlier_id <- Gilliamella_apicola_vs_recodes_only$id1
Gilliamella_apicola_vs_recodes_only$outlier_id[which(Gilliamella_apicola_vs_recodes_only$id2 %in% apicola_recodes)] <- Gilliamella_apicola_vs_recodes_only[which(Gilliamella_apicola_vs_recodes_only$id2 %in% apicola_recodes), 'id2']
aggregate(ani ~ outlier_id, data = Gilliamella_apicola_vs_recodes_only, FUN = max)


# Bartonella
within_Bartonella <- subset_fastani(fastani_out, 'Bartonella_apis')
round(mean(within_Bartonella$ani), 1)
round(sd(within_Bartonella$ani), 1)


# Bifidobacterium_asteroides outliers
within_asteroides <- subset_fastani(fastani_out, 'Bifidobacterium_asteroides')
asteroides_outliers <- c('GCA_003202695.1')
within_asteroides_outliers_id1_i <- which(within_asteroides$id1 %in% asteroides_outliers)
within_asteroides_outliers_id2_i <- which(within_asteroides$id2 %in% asteroides_outliers)
within_asteroides_outliers_i <- c(setdiff(within_asteroides_outliers_id2_i, within_asteroides_outliers_id1_i),
                                  setdiff(within_asteroides_outliers_id1_i, within_asteroides_outliers_id2_i))
asteroides_vs_outliers_only <- within_asteroides[within_asteroides_outliers_i, ]
asteroides_vs_outliers_only$outlier_id <- asteroides_vs_outliers_only$id1
asteroides_vs_outliers_only$outlier_id[which(asteroides_vs_outliers_only$id2 %in% asteroides_outliers)] <- asteroides_vs_outliers_only[which(asteroides_vs_outliers_only$id2 %in% asteroides_outliers), 'id2']

within_asteroides_max_raw1 <- aggregate(ani ~ id1, data = within_asteroides, FUN = max)
within_asteroides_max_raw2 <- aggregate(ani ~ id2, data = within_asteroides, FUN = max)
colnames(within_asteroides_max_raw1) <- c('id', 'ani')
colnames(within_asteroides_max_raw2) <- c('id', 'ani')

within_asteroides_max <- aggregate(ani ~ id, data = rbind(within_asteroides_max_raw1, within_asteroides_max_raw2), FUN = max)
