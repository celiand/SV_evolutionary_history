safe_read <- function(f) {
  if (!file.exists(f) || file.info(f)$size == 0) {
    message(f, " is missing or empty")
    return(NULL)   # or return an empty data.frame()
  }
  read.table(f, header = FALSE, stringsAsFactors = FALSE)
}

g1 <- safe_read("group1_ids.txt")
g2 <- safe_read("group2_ids.txt")

if (is.null(g1) || is.null(g2)) {
  ## cant run the analysis
}else{
  
  
  library(data.table)
  library(ggplot2)
  
  # Load sample IDs
  id <- fread("region_ibs.mdist.id", header = FALSE, col.names = c("FID","IID"))
  M  <- as.matrix(fread("region_ibs.mdist", header = FALSE))
  stopifnot(nrow(M) == nrow(id), ncol(M) == nrow(id))
  rownames(M) <- id$IID
  colnames(M) <- id$IID
  
  # Define groups (provide two text files with one sample ID per line)
  g1_ids <- readLines("group1_ids.txt")  # Individuals carrying AA genotype
  g2_ids <- readLines("group2_ids.txt")  # Individuals carrying aa genotype
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
  
  
  ### save data to plot
  filenameibdfile<-paste("/mnt/SCRATCH/cedi/phDSalmon/slimfolder/vcfkeptresults/allpop",filename,"_idbdata.txt",sep="")
  
  write.table(pairs_dt,file=filenameibdfile,row.names = FALSE)
  
}