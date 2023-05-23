rm(list = ls(all.names = TRUE))

# Dataset samples
Ellegaard.2019_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_ids.txt", stringsAsFactors = FALSE)$V1
Ellegaard.2020_samples <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_ids.txt", stringsAsFactors = FALSE)$V1

sample_names_2019 <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2019_SRR_to_sample.txt",
                                header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sample_names_2019) <- c("SRR", "sample")
rownames(sample_names_2019) <- sample_names_2019$SRR

sample_names_2019$colony <- gsub("Y.*$", "", sample_names_2019$sample)

sample_names_2019$location <- gsub("_.*$", "", sample_names_2019$sample)
sample_names_2019$location <- gsub("^.*Y", "", sample_names_2019$location)

sample_names_2019$age_raw <- gsub(".*_", "", sample_names_2019$sample)
sample_names_2019$age <- NA
sample_names_2019[grep("W", sample_names_2019$age_raw), "age"] <- "Old"
sample_names_2019[grep("F", sample_names_2019$age_raw), "age"] <- "Middle-aged"
sample_names_2019[grep("N", sample_names_2019$age_raw), "age"] <- "Young"


sample_names_2020 <- read.table("/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard2020_SRR_to_sample.txt",
                                header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(sample_names_2020) <- c("SRR", "sample")
rownames(sample_names_2020) <- sample_names_2020$SRR

sample_names_2020$Apiary <- NA
sample_names_2020$Apiary[grep("Ai", sample_names_2020$sample)] <- "Ai"
sample_names_2020$Apiary[grep("Iu", sample_names_2020$sample)] <- "Iu"



# Fill in combined metadata table.
sample_metadata <- data.frame(matrix(NA, nrow = length(c(Ellegaard.2019_samples, Ellegaard.2020_samples)), ncol = 5))
rownames(sample_metadata) <- c(Ellegaard.2019_samples, Ellegaard.2020_samples)
colnames(sample_metadata) <- c("Name", "Country", "Apiary", "Year", "Age")

# Name
sample_metadata[Ellegaard.2019_samples, "Name"] <- sample_names_2019[Ellegaard.2019_samples, "sample"]
sample_metadata[Ellegaard.2020_samples, "Name"] <- sample_names_2020[Ellegaard.2020_samples, "sample"]


# Country
sample_metadata[Ellegaard.2019_samples, "Country"] <- "Switzerland"
sample_metadata[Ellegaard.2020_samples, "Country"] <- "Japan"


# Year
sample_metadata[sample_names_2019[which(sample_names_2019$location == "1"), "SRR"], "Year"] <- "2019_Y1"
sample_metadata[sample_names_2019[which(sample_names_2019$location == "2"), "SRR"], "Year"] <- "2019_Y2"
sample_metadata[Ellegaard.2020_samples, "Year"] <- "2020"


# Apiary
sample_metadata[sample_names_2019[which(sample_names_2019$colony == "Dr"), "SRR"], "Apiary"] <- "Les Droites"
sample_metadata[sample_names_2019[which(sample_names_2019$colony == "Gr"), "SRR"], "Apiary"] <- "Grammont"
sample_metadata[sample_names_2020[which(sample_names_2020$Apiary == "Ai"), "SRR"], "Apiary"] <- "Ai"
sample_metadata[sample_names_2020[which(sample_names_2020$Apiary == "Iu"), "SRR"], "Apiary"] <- "Iu"


# Age
sample_metadata[sample_names_2019$SRR, "Age"] <- sample_names_2019$age
sample_metadata[Ellegaard.2020_samples, "Age"] <- "Young"


write.table(x = sample_metadata, file = "/data1/gdouglas/projects/honey_bee/combined_Ellegaard.2019.2020/mapfiles/Ellegaard.2019.2020_metadata.tsv",
            row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
