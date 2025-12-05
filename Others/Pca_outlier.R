### Make a genome wide PCA using Bayestyper_SVs
library(tidyverse)

#### first load samples name and populations (can be extracted from vcf)

AM_SAMP_D<-read.table("Data/Farmed_American_sample_SNP_Compat_2.txt")
AM_SAMP_D$pop<-"Domesticated_American"

EU_SAMP_D<-read.table("Data/Farmed_European_sample_SNP_Compat_2.txt")
EU_SAMP_D$pop<-"Domesticated_European"

AM_SAMP_W<-read.table("Data/Wild_American_sample_SNP_Compat_2.txt")
AM_SAMP_W$pop<-"Wild_American"

EU_SAMP_W<-read.table("Data/Wild_European_sample_SNP_Compat.txt")
EU_SAMP_W$pop<-"Wild_European"

popdata<-rbind.data.frame(AM_SAMP_D,EU_SAMP_D,AM_SAMP_W,EU_SAMP_W)
colnames(popdata)<-c("IID","pop")

##then get pca values (from plink2 --pca command)
region_pca <-data.table::fread("Data/pca_SV_res.eigenvec")  #for entire PCA plot

region_pca_var <-data.table::fread("Data/pca_SV_res.eigenval") 


colnames(region_pca_var)<-"eigenval"
region_pca_var<- region_pca_var %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 

prop2<-as.numeric(as.character(region_pca_var[2,2]))*100
prop1<-as.numeric(as.character(region_pca_var[1,2]))*100

#input the populations file
pops_file <- popdata


#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged <- merge(region_pca, pops_file, by.x = "#IID",by.y = "IID")

#Plot

ylabel<-paste("PC2 - ",prop2,"%",sep="")
xlabel<-paste("PC1 - ",prop1,"%",sep="")

ggplot(PCA_file_merged,aes(x = PC1, y = PC2, group=pop)) + geom_point((aes(colour = pop)), size = 4, show.legend = TRUE) + 
  theme_bw(25) + labs(x = xlabel, y = ylabel) +scale_color_manual(values=c("#33a02c","#1f78b4","#b2df8a","#a6cee3"))