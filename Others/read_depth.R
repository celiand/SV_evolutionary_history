library(tidyverse)


output_file <- "rrd_matrix.tsv"


base_dir <- "Data/mosdepth_results"
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

rrd_list <- list()
chrom_map <- tibble(
  chrom = paste0("NC_059", sprintf("%03d", 442:470), ".1"),
  ssa = paste0("ssa", sprintf("%02d", 1:29))
)

##process samples
for (sample_dir in sample_dirs) {
  sample <- basename(sample_dir)
  message("Processing: ", sample)
  
  region_files <- list.files(sample_dir, pattern = "ssa.*\\.regions\\.bed\\.gz$", full.names = TRUE)
  if (length(region_files) == 0) next
  
  df <- map_dfr(region_files, ~{
    read_tsv(.x, col_names = c("chrom", "start", "end", "depth"), show_col_types = FALSE) %>%
      left_join(chrom_map, by = "chrom") %>%
      mutate(chrom = coalesce(ssa, chrom)) %>%
      select(-ssa)
  })
  
  genome_mean <- mean(df$depth, na.rm = TRUE)
  
  df <- df %>%
    mutate(sample = sample,
           rrd = depth / genome_mean) %>%
    select(chrom, start, end, sample, rrd)
  
  rrd_list[[sample]] <- df
}


## save data
rrd_long <- bind_rows(rrd_list) %>%
  arrange(chrom, start)


write_tsv(rrd_long, "rrd_long.tsv")


library(tidyverse)


## we can also make a matrix
rrd_matrix <- rrd_long %>%
  group_by(chrom) %>%
  group_split() %>%
  map(~ pivot_wider(.x, names_from = sample, values_from = rrd, values_fill = 0)) %>%
  bind_rows() %>%
  arrange(chrom, start)

rrd_matrix_rounded <- rrd_matrix %>%
  mutate(across(-c(chrom, start, end), ~ round(.x, 2)))


rrd_matrix_rounded %>%
  group_by(chrom) %>%
  group_split() %>%
  walk(~ write_tsv(.x, paste0("rrd_matrix_", unique(.x$chrom), ".tsv")))
