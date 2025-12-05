## ---- First format data to make overlap ----

## input pca region
regionfile<-read.table(file="Data/interesting_region_bis.txt",header=TRUE)

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chrom <- sprintf("ssa%02d", 1:29)

##define range to look at around Sv boundary
range<-100000

### generate random boundaries: random location in the genome that we will compare to actual SV boundaries

nbfakeboundary<-1000
allfakeboundary<-c()
for (i in 1:nbfakeboundary){
  i <- sample(1:29, 1)
  chr <- chrom[i]
  chr_len <- lengthvector[i]
  
  # choose a valid random point that respects the range constraints
  pos <- sample((range + 1):(chr_len - range - 1), 1)
  
  lower <- pos - range
  upper <- pos + range
  
  result <- data.frame(chromosome = chr, start = lower, end = upper)
  allfakeboundary<-rbind.data.frame(allfakeboundary,result,stringsAsFactors = FALSE)
}

allfakeboundary$status<-"ramdom"

##now make a same format table for real SV
empiricalboundary_start<-data.frame(chromosome=regionfile$CHROM,start=regionfile$Manual_start-range,end=regionfile$Manual_start+range)
empiricalboundary_start$status<-"notrecurrent"
empiricalboundary_start[c(6,15,17,22),"status"]<-"recurrent" ### change the status of recurrent Sv , from the Ibs analysis


empiricalboundary_end<-data.frame(chromosome=regionfile$CHROM,start=regionfile$Manual_end-range,end=regionfile$Manual_end+range)
empiricalboundary_end$status<-"notrecurrent"
empiricalboundary_end[c(6,15,17,22),"status"]<-"recurrent" ### change the status of recurrent Sv , from the Ibs analysis

allboundaries<-rbind.data.frame(allfakeboundary,empiricalboundary_start,empiricalboundary_end,stringsAsFactors = FALSE)


### export data to make an overlap using bedtools
options("scipen"=100, "digits"=4) ##disable scientific notation for later exportation

write.table(allboundaries,file="allboundariSV_testTE.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")


## ---- Then load overlap data for analysis ----

overlapboundaryTE<-read.table(file="Data/allboundariSV_testTE_overlap_repeat.txt",header=FALSE,sep="\t")

colnames(overlapboundaryTE)<-c("chrboundary","startboundary","endboundary","status","chromTE","startTE","endTE","typeTE","divergence")
overlapboundaryTE$IDTE<-paste(overlapboundaryTE$chromTE,overlapboundaryTE$startTE,overlapboundaryTE$endTE,overlapboundaryTE$typeTE,sep="_")
overlapboundaryTE<-overlapboundaryTE[overlapboundaryTE$typeTE!="Segmental duplication",] #### note: here we can change  the != to == to have only segmental duplications

library(tidyverse)

boundary_summary <- overlapboundaryTE %>%
  mutate(IDboundary = paste(chrboundary, startboundary, endboundary, sep = "_")) %>%
  group_by(IDboundary, status) %>%
  summarise(nbofTE = n_distinct(IDTE), .groups = "drop")%>%
  select(IDboundary, nbofTE, status)  

# number of expected boundary per status
required_counts <- c(ramdom = 1000, notrecurrent = 44, recurrent = 8)

# Count how many we actually have
current_counts <- table(boundary_summary$status)

rows_to_add <- lapply(names(required_counts), function(status) {
  current <- ifelse(status %in% names(current_counts), current_counts[status], 0)
  missing <- required_counts[status] - current
  
  if (missing > 0) {
    tibble(
      IDboundary = paste0("leftover_", seq_len(missing), "_", status),
      nbofTE = 0,
      status = status
    )
  } else {
    NULL
  }
}) %>% bind_rows()

boundary_final <- bind_rows(boundary_summary, rows_to_add)

boundary_final <- boundary_final %>%
  mutate(status = recode(status, ramdom = "random")) %>%
  mutate(status = factor(status, levels = c("recurrent", "notrecurrent", "random")))


### ---- make a combined plot ---- 

df_box <- boundary_final %>% filter(status %in% c("recurrent", "notrecurrent"))
df_density <- boundary_final %>% filter(status == "random")


library(ggplot2)
library(patchwork)
library(cowplot)

mycols <- c("random" = "#4477AA", "notrecurrent" = "#EE6677", "recurrent" = "#228833")

p_box <- ggplot(df_box, aes(x = nbofTE, y = status, fill = status)) +
  scale_x_continuous(limits = c(0, 2000))+
  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = 21) +
  scale_fill_manual(values = mycols) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  ylab("") +
  xlab("Number of TEs") +
  ggtitle("Recurrent & Non-recurrent")



p_density <- ggplot(df_density, aes(x = nbofTE, fill = status)) +
  geom_density(alpha = 0.4, color = "black") +
  scale_fill_manual(values = mycols) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  ylab("Density") +
  xlab("") +
  ggtitle("Random boundaries")+scale_x_continuous(limits = c(0, 2000))

combined_plot <- plot_grid(p_box,p_density  , ncol = 1, align = "v")
combined_plot


### ---- find statistical groups ---

library(FSA)
library(multcompView)
library(dplyr)

# Kruskal–Wallis
kw <- kruskal.test(nbofTE ~ status, data = boundary_final)

# Dunn post-hoc
dunn_res <- dunnTest(nbofTE ~ status, data = boundary_final, method = "bh")
dunn
# Convert p-values into pairwise compact letters
pvals <- dunn_res$res %>%
  select(Comparison, P.adj)


# Extract p-values into a named vector
pvec <- pvals$P.adj
names(pvec) <- gsub(" ", "", pvals$Comparison)   # remove spaces => "A-B"
pvec

letters <- multcompLetters(pvec)$Letters
letters