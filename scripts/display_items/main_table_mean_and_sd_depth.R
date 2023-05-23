rm(list = ls(all.names = TRUE))

# Read in SRRs per dataset.
datasets <- c('Ellegaard2019', 'Ellegaard2020', 'Sun2022', 'Wu2021', 'Zhang2022')
SRRs <- list()
for (dataset in datasets) {
  SRRs[[dataset]] <- read.table(paste('/data1/gdouglas/projects/bee_microbiome_figshare/datasets/SRRs/', dataset, '_SRRs.txt.gz', sep = ''),
                                stringsAsFactors = FALSE)$V1
}

# Read in number of reads per FASTQ.
num_reads <- read.delim('/data1/gdouglas/projects/bee_microbiome_figshare/datasets/fastq_num_reads.txt.gz',
                        header = FALSE, stringsAsFactors = FALSE, sep = ' ', row.names = 1)
colnames(num_reads) <- 'count'


# Run quick sanity check that all paired-end FASTQs have the same number of forward and reverse reads.
paired_R1_fastqs <- grep("_paired.1.fastq.gz", rownames(num_reads), value = TRUE)
paired_R2_fastqs <- gsub('_paired.1.fastq.gz', '_paired.2.fastq.gz', paired_R1_fastqs)
identical(num_reads[paired_R1_fastqs, 'count'], num_reads[paired_R2_fastqs, 'count'])
# TRUE


# Sum total read counts per SRR, and then compute mean and SD per sample.
num_reads$SRR <- gsub('_.*$', '', rownames(num_reads))

num_reads_sum <- aggregate(count ~ SRR, data = num_reads, FUN = sum)

num_reads_sum$study <- NA
rownames(num_reads_sum) <- num_reads_sum$SRR
for (dataset in datasets) {
  num_reads_sum[SRRs[[dataset]], 'study'] <- dataset
}

num_reads_mean <- aggregate(count ~ study, data = num_reads_sum, FUN = mean)
num_reads_sd <- aggregate(count ~ study, data = num_reads_sum, FUN = sd)
rownames(num_reads_mean) <- num_reads_mean$study
rownames(num_reads_sd) <- num_reads_sd$study

summary_tab <- data.frame(datasets = names(SRRs),
                          num_samples = sapply(SRRs, length),
                          mean_num_millions = round(num_reads_mean[names(SRRs), 'count'] / 1e6, 1),
                          sd_num_millions = round(num_reads_sd[names(SRRs), 'count'] / 1e6, 1))

# Also add in mean and sd as text in clean format.
summary_tab$mean_and_sd_char <- paste(as.character(summary_tab$mean_num_millions),
                                      ' (',
                                      as.character(summary_tab$sd_num_millions),
                                      ')',
                                      sep = '')

write.table(x = summary_tab,
            file = '/home/gdouglas/scripts/bee_microbome_allelic_div/display_items/Table1_raw.tsv',
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
