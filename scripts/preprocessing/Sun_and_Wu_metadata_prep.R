# Clean-up and Sun and Wu dataset metadata files and write out clean tables.
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
