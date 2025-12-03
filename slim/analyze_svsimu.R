### the aim is to have an overview of the retention of SV in the different simulations types

library(tidyverse)
library(dplyr)
simu_table<-read.table(file="all_vcf_info.txt",header=TRUE,fill=TRUE)
simu_table_clean<-simu_table[!is.na(simu_table$Run),]

### keep only the SV
simu_table_clean_SV<-simu_table_clean[simu_table_clean$Mutation_type==2,]

siminfopop<-simu_table_clean_SV %>% group_by(Type_simu,Run) %>% summarize(count_frequency_gt_0 = sum(Frequency > 0), .groups = "drop")


library(ggplot2)


simupop_retain<-siminfopop %>% group_by(count_frequency_gt_0,Type_simu) %>%
  summarize(nb_ofocc = n())


legend_order <- c("4", "3", "2", "1")

simupop_retain<-simupop_retain[simupop_retain$count_frequency_gt_0%in%c("1","2","3","4"),]



simupop_retain_proportion <- simupop_retain %>%
  group_by(Type_simu) %>%
  mutate(proportion = nb_ofocc / sum(nb_ofocc))

counts <- simupop_retain %>%
  group_by(Type_simu) %>%
  summarise(total_count = sum(nb_ofocc))


simupop_retain_proportion$count_frequency_gt_0 <- factor(
  simupop_retain_proportion$count_frequency_gt_0,
  levels = c(4, 3, 2, 1) # Specify the desired order
)


### for pop = 2 split same continent separately

simu_table_clean_SV$Population<-as.numeric(as.character(simu_table_clean_SV$Population))

siminfopopbis <- simu_table_clean_SV %>%
  group_by(Type_simu, Run) %>%
  summarize(
    count_frequency_gt_0 = sum(Frequency > 0),
    sumpop = sum(Population[Frequency > 0]),
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      count_frequency_gt_0 == 2 & sumpop %in% c(6, 8) ~ "Same_continent",
      TRUE ~ "Diff_continent"
    )
  )



# Merge label info back in
simupop_retain_labeled <- siminfopopbis %>%
  group_by(Type_simu, count_frequency_gt_0, label) %>%
  summarize(nb_ofocc = n(), .groups = "drop")

# Compute proportions
simupop_retain_proportion_bis <- simupop_retain_labeled %>%
  group_by(Type_simu) %>%
  mutate(proportion = nb_ofocc / sum(nb_ofocc))

# Create a new fill variable
simupop_retain_proportion_bis <- simupop_retain_proportion_bis %>%
  mutate(
    fill_category = case_when(
      count_frequency_gt_0 == 2 ~ paste0("2_", label),
      TRUE ~ as.character(count_frequency_gt_0)
    ),
    count_frequency_gt_0 = factor(count_frequency_gt_0, levels = c(4, 3, 2, 1))
  )



library(ggpattern)
library(ggplot2)




#### ok now distinguish SV fixed in all 4 populations vs Sv not fixed in all 4 populations (polymorphic) 

siminfopoptre <- simu_table_clean_SV %>%
  group_by(Type_simu, Run) %>%
  summarize(
    count_frequency_gt_0 = sum(Frequency > 0),
    sumpop = sum(Population[Frequency > 0]),
    all_freq_equal_1 = all(Frequency[Frequency > 0] == 1),
    .groups = "drop"
  ) %>%
  mutate(
    label = case_when(
      count_frequency_gt_0 == 2 & sumpop %in% c(6, 8) ~ "Same_continent",
      TRUE ~ "Diff_continent"
    )
  )

#### we are only interest in cases where there are 4 pop
siminfopoptre[siminfopoptre$count_frequency_gt_0!=4,"all_freq_equal_1"]<-TRUE


simupop_retain_labeled <- siminfopoptre %>%
  group_by(Type_simu, count_frequency_gt_0, label, all_freq_equal_1) %>%
  summarize(nb_ofocc = n(), .groups = "drop")

# Compute proportions
simupop_retain_proportion_tre <- simupop_retain_labeled %>%
  group_by(Type_simu) %>%
  mutate(proportion = nb_ofocc / sum(nb_ofocc))


simupop_retain_proportion_tre <- simupop_retain_proportion_tre %>%
  mutate(
    fill_category = case_when(
      count_frequency_gt_0 == 2 ~ paste0("2_", label),
      count_frequency_gt_0 == 4 & all_freq_equal_1 ~ "4_fixed",
      count_frequency_gt_0 == 4 & !all_freq_equal_1 ~ "4_polymorphic",
      TRUE ~ as.character(count_frequency_gt_0)
    ),
    count_frequency_gt_0 = factor(count_frequency_gt_0, levels = c(4, 3, 2, 1))
  )

ggplot(simupop_retain_proportion_tre, aes(
  x = Type_simu,
  y = proportion,
  fill = fill_category,
  pattern = fill_category
)) +
  geom_bar_pattern(
    stat = "identity",
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.3,
    pattern_angle = 45,
    pattern_fill = "black",
    pattern_colour = "black"
  )+scale_fill_manual(
    values = c(
      "4_fixed" = "#006d2c",
      "4_polymorphic" = "#006d2c",
      "3" = "#238b45",
      "2_Same_continent" = "#66c2a4",
      "2_Diff_continent" = "#66c2a4",
      "1" = "#b2e2e2"
    )
  ) +
  scale_pattern_manual(
    values = c(
      "4_fixed" = "none",
      "4_polymorphic" = "stripe",
      "3" = "none",
      "2_Same_continent" = "circle",
      "2_Diff_continent" = "none",
      "1" = "none"
    )
  )+ theme_bw(20) +theme(axis.text.x = element_text(angle = 45, hjust=1))



## only 1.1 scenarios, and keep only 1000 rec rate 

simupop_retain_proportion_tre_filter <- simupop_retain_proportion_tre %>% filter(str_starts(Type_simu, "1.1")) %>% filter(!str_detect(Type_simu, "(^|_)500(_|$)|(^|_)1500(_|$)"))

ggplot(simupop_retain_proportion_tre_filter, aes(
  x = Type_simu,
  y = proportion,
  fill = fill_category,
  pattern = fill_category
)) +
  geom_bar_pattern(
    stat = "identity",
    pattern_density = 0.2,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.3,
    pattern_angle = 45,
    pattern_fill = "black",
    pattern_colour = "black"
  )+scale_fill_manual(
    values = c(
      "4_fixed" = "#006d2c",
      "4_polymorphic" = "#006d2c",
      "3" = "#238b45",
      "2_Same_continent" = "#66c2a4",
      "2_Diff_continent" = "#66c2a4",
      "1" = "#b2e2e2"
    )
  ) +
  scale_pattern_manual(
    values = c(
      "4_fixed" = "none",
      "4_polymorphic" = "stripe",
      "3" = "none",
      "2_Same_continent" = "circle",
      "2_Diff_continent" = "none",
      "1" = "none"
    )
  )+ theme_bw(20) +theme(axis.text.x = element_text(angle = 45, hjust=1))






##### make a random samples of 12 vcf for each cases
keptsimu<-unique(simupop_retain_proportion_tre_filter$Type_simu)

firstsubset<-siminfopoptre[siminfopoptre$count_frequency_gt_0==4 & siminfopoptre$Type_simu%in%keptsimu & siminfopoptre$all_freq_equal_1==FALSE,]

sampledkept <- NULL
for( i in keptsimu){
  subsettypesimu<-firstsubset[firstsubset$Type_simu==i,]
  if(nrow(subsettypesimu)==0){
    
  }else if(nrow(subsettypesimu)==12){
    sampledkept <- rbind(sampledkept, subsettypesimu)
  }else{
    sampledkept <- rbind(sampledkept, subsettypesimu[sample(nrow(subsettypesimu), 12), ])
  }
}

write.table(sampledkept,file="sampledkept.txt",row.names = FALSE)

head(sampledkept)