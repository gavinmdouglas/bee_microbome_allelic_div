# Clean-up dataset metadata files and write out clean tables.
# The key action is to link the sample metadata to each SRR id.

# For Sun dataset, first read in Supp Table with sample metadata.
# Note that this was restricted to A. mellifera samples only.
# For this dataset, there actually isn't any variation across the provided 
# metadata, but cleaning up just to make this clear for the future.
rm(list = ls(all.names = TRUE))

Sun_supp <- read.table('/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/mapfiles/raw_metadata/Sun_raw_metadata/40168_2022_1268_MOESM2_ESM.txt',
                       header = TRUE, sep = '\t', stringsAsFactors = FALSE)
Sun_supp <- Sun_supp[which(Sun_supp$Sample_ID != ''), ]

# Then read in table from ENA ('filereport') with SRR ids and sample ids.
Sun_filereport <- read.table('/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/mapfiles/raw_metadata/Sun_raw_metadata/filereport_read_run_PRJNA787435_tsv.txt',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)
rownames(Sun_filereport) <- Sun_filereport$sample_alias

Sun_supp_orig_colnames <- colnames(Sun_supp)
Sun_supp$SRR <- Sun_filereport[Sun_supp$Sample_ID, 'run_accession']
Sun_supp <- Sun_supp[, c('SRR', Sun_supp_orig_colnames)]

write.table(x = Sun_supp,
            file = '/data1/gdouglas/projects/bee_microbiome_figshare/datasets/metadata/Sun2022_metadata.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


# Then for Wu dataset - simplify filereport table to just be SRRs, sample ids, and sample groups.
rm(list = ls(all.names = TRUE))

Wu_filereport <- read.table('/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/mapfiles/raw_metadata/Wu_raw_metadata/filereport_read_run_PRJNA645015_tsv.txt',
                            header = TRUE, sep = '\t', stringsAsFactors = FALSE)

Wu_metadata <- Wu_filereport[, c('run_accession', 'sample_alias')]

colnames(Wu_metadata) <- c('SRR', 'Sample_ID')

Wu_metadata$Group <- gsub('\\d+$', '', Wu_metadata$Sample_ID)

write.table(x = Wu_metadata,
            file = '/data1/gdouglas/projects/bee_microbiome_figshare/datasets/metadata/Wu2021_metadata.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


# For Zhang dataset, just re-arrange columns that had previously been created in excel by hand.
rm(list = ls(all.names = TRUE))
Zhang_metadata <- read.table('/data1/gdouglas/projects/honey_bee/Chinese_and_Ellegaard/mapfiles/raw_metadata/Zhang_raw_metadata/Zhang2022_learning_PRJNA670603_metadata.txt',
                             header = TRUE, sep = '\t', stringsAsFactors = FALSE)

Zhang_metadata <- Zhang_metadata[, c('SRR', 'Sample', 'Day', 'Treatment')]
colnames(Zhang_metadata)[2] <- 'Sample_ID'

write.table(x = Zhang_metadata,
            file = '/data1/gdouglas/projects/bee_microbiome_figshare/datasets/metadata/Zhang2022_metadata.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


# Already had metadata prepped for each Ellegaard dataset, so just split these into two separate tables (just for consistency).
rm(list = ls(all.names = TRUE))
Ellegaard_metadata <- read.table('/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard.2019.2020_metadata.tsv',
                                 header = TRUE, sep = '\t', stringsAsFactors = FALSE)

colnames(Ellegaard_metadata) <- c('SRR', 'Sample_ID', 'Country', 'Apiary', 'Year', 'Age')

Ellegaard_metadata_2019 <- Ellegaard_metadata[which(Ellegaard_metadata$Country == 'Switzerland'), ]
Ellegaard_metadata_2020 <- Ellegaard_metadata[which(Ellegaard_metadata$Country == 'Japan'), ]

write.table(x = Ellegaard_metadata_2019,
            file = '/data1/gdouglas/projects/bee_microbiome_figshare/datasets/metadata/Ellegaard2019_metadata.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')

write.table(x = Ellegaard_metadata_2020,
            file = '/data1/gdouglas/projects/bee_microbiome_figshare/datasets/metadata/Ellegaard2020_metadata.tsv',
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
