genotypeofSV<-read.table(file="Data/bayestyper_ALL_SV_filtered_renamed_sorted.recode.vcf.gz_258.stat")
samplefile<-read.table(file="Data/samplebttest.txt") ## file with sample named

new_colnames <- c("CHROM", samplefile[[1]])
colnames(genotypeofSV) <- new_colnames



library(stringr)

path <- "Data/ibsfolder_bt"
files <- list.files(path, full.names = FALSE)

# Extract the base region names
regions_mdist     <- str_match(files, "^(.*)_rg_ibs\\.mdist$")[,2]
regions_mdist_id  <- str_match(files, "^(.*)_rg_ibs\\.mdist\\.id$")[,2]

# Keep only unique regions that have BOTH .mdist and .mdist.id
unique_names <- intersect(na.omit(regions_mdist), na.omit(regions_mdist_id))



library(data.table)

process_region <- function(region, mdist_path, genotype_table) {
  # Build filenames
  file_id   <- file.path(mdist_path, paste0(region, "_rg_ibs.mdist.id"))
  file_mdist <- file.path(mdist_path, paste0(region, "_rg_ibs.mdist"))
  
  # Read IBS files
  id <- fread(file_id, header = FALSE, col.names = c("FID","IID"))
  M  <- as.matrix(fread(file_mdist, header = FALSE))
  stopifnot(nrow(M) == nrow(id), ncol(M) == nrow(id))
  rownames(M) <- id$IID
  colnames(M) <- id$IID
  
  # --- Get SV position from region name ---
  pos <- as.integer(strsplit(region, "_")[[1]][2])  # e.g. "ssa17_7147718_..." → 7147718
  
  # Extract genotype row from genotypeofSV
  sv_row <- genotype_table[genotype_table$CHROM == pos, ]
  if (nrow(sv_row) == 0) stop(paste("No genotype found for", region))
  
  # Melt into long format: IID + genotype
  geno_df <- data.frame(IID = colnames(sv_row)[-1],
                        GT = as.character(t(sv_row[,-1, drop=FALSE])))
  
  # Define groups (0/0 vs 1/1)
  g1 <- geno_df$IID[geno_df$GT == "0/0"]
  g2 <- geno_df$IID[geno_df$GT == "1/1"]
  
  # Intersect with available IDs
  g1 <- intersect(g1, id$IID)
  g2 <- intersect(g2, id$IID)
  
  # Convert upper triangle of distance matrix into long format
  upper_idx <- which(upper.tri(M), arr.ind = TRUE)
  pairs_dt <- data.table(
    i = rownames(M)[upper_idx[,1]],
    j = colnames(M)[upper_idx[,2]],
    dist = M[upper_idx]
  )
  
  # Assign pair types
  pairs_dt[, pair_type := fifelse(i %in% g1 & j %in% g1, "G1-G1",
                                  fifelse(i %in% g2 & j %in% g2, "G2-G2",
                                          fifelse((i %in% g1 & j %in% g2) | (i %in% g2 & j %in% g1), "G1-G2", NA_character_)))]
  pairs_dt <- pairs_dt[!is.na(pair_type)]
  
  # Keep only regions that have all 3 groups
  if (length(unique(pairs_dt$pair_type)) < 3) {
    return(NULL)  # skip this region
  }
  
  pairs_dt[, region := region]
  
  return(pairs_dt)
}



mdist_path <- "Data/ibsfolder_bt"

# Example for one region
region <- unique_names[1]
pairs_dt <- process_region(region, mdist_path, genotypeofSV)

# Or loop over all
all_pairs <- rbindlist(lapply(unique_names, function(r) process_region(r, mdist_path, genotypeofSV)))


write.table(all_pairs,"allbt_ibsdata.txt",row.names=FALSE)

dataall<-read.table(file="allbt_ibsdata.txt",header=TRUE)


library(ggplot2)
library(dplyr)


### filter region with too low number of SNPs (less range of values)
all_pairs_filtered <- dataall %>%
  group_by(region) %>%
  filter(n_distinct(dist) >= 30) %>%
  ungroup()

unique(all_pairs_filtered$region)


###plot one region
pairs_dt<-dataall[dataall$region=="ssa29_25042963_25050093",]


ggplot(pairs_dt, aes(x = dist, fill = pair_type)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 40) +
  labs(x = "Distance (1 - IBS)", y = "Number of pairs", fill = "Pair type") +
  theme_classic()