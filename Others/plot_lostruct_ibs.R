library(data.table)
library(ggplot2)
# Load sample IDs

path <- "Data/ibsfolder/"

# List all files ending with _rg.bim
files <- list.files(path, pattern = "_rg\\.bim$", full.names = FALSE)

# Remove the trailing "_rg.bim"
names_clean <- sub("_rg\\.bim$", "", files)

# Keep unique values
unique_names <- unique(names_clean)

unique_names

### plot a single region
### look recurrent : names 6, 15, 17, 22

region<-unique_names[10]
filneameibd<-paste("Data/ibsfolder/",region,"_rg_ibs.mdist.id",sep="")
filneamemist<-paste("Data/ibsfolder/",region,"_rg_ibs.mdist",sep="")
id <- fread(filneameibd, header = FALSE, col.names = c("FID","IID"))
M  <- as.matrix(fread(filneamemist, header = FALSE))
stopifnot(nrow(M) == nrow(id), ncol(M) == nrow(id))
rownames(M) <- id$IID
colnames(M) <- id$IID

filenamepcagt<-paste("Data/sprgi_gt_pca_files/",region,".txt",sep="")
region_data<-read.table(file=filenamepcagt,header=TRUE)



# Define groups (provide two text files with one sample ID per line)
g1_ids <- region_data[region_data$GT=="I/I","IID"]  # Individuals carrying AA genotype
g2_ids <- region_data[region_data$GT=="NI/NI","IID"]  # Individuals carrying aa genotype
g1 <- intersect(g1_ids, id$IID)
g2 <- intersect(g2_ids, id$IID)
# Convert upper triangle of distance matrix into long format
upper_idx <- which(upper.tri(M), arr.ind = TRUE)
pairs_dt <- data.table(
  i = rownames(M)[upper_idx[,1]],
  j = colnames(M)[upper_idx[,2]],
  dist = M[upper_idx]
)
# Assign pair type (G1-G1, G2-G2, or G1-G2)
pairs_dt[, pair_type := fifelse(i %in% g1 & j %in% g1, "G1-G1",
                                fifelse(i %in% g2 & j %in% g2, "G2-G2",
                                        fifelse((i %in% g1 & j %in% g2) | (i %in% g2 & j %in% g1), "G1-G2", NA_character_)))]
pairs_dt <- pairs_dt[!is.na(pair_type)]
# Plot histograms
ggplot(pairs_dt, aes(x = dist, fill = pair_type)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 40) +
  labs(x = "Distance (1 - IBS)", y = "Number of pairs", fill = "Pair type") +
  theme_classic()
