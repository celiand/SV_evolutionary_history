library(scales)
library(ggplot2)

## ---- plot size distribution of manta SV and bayestyper SV  ----

FE_SV_manta<-read.table(file="Data/SV_info_table_Farmed_Europe.txt",header=TRUE)

FE_SV_manta$size<-FE_SV_manta$END-FE_SV_manta$POS

FE_SV_bayestyper<-read.table(file="Data/bayestyper_SV_info_table_Farmed_European.txt",header=TRUE)

FE_SV_bayestyper$size<-FE_SV_bayestyper$END-FE_SV_bayestyper$POS

summary(FE_SV_manta[FE_SV_manta$size>50,"size"])

summary(FE_SV_bayestyper[FE_SV_bayestyper$size>50,"size"])


ggplot(FE_SV_bayestyper[FE_SV_bayestyper$size>50,], aes(x = size)) + 
  geom_histogram(color = "black", fill = "dodgerblue3", alpha = 0.7, bins = 30) + 
  scale_x_log10(labels = comma,limits = c(50, 100000))+theme_bw()


ggplot(FE_SV_manta[FE_SV_manta$size>50,], aes(x = size)) + 
  geom_histogram(color = "black", fill = "dodgerblue3", alpha = 0.7, bins = 30) + 
  scale_x_log10(labels = comma,limits = c(50, 100000))+theme_bw()


## ---- extract location of size peak SV  ----


putative_tc1<-FE_SV_bayestyper[FE_SV_bayestyper$size>=1432 & FE_SV_bayestyper$size<=1436,]
options("scipen"=100, "digits"=4) ##disable scientific notation
putative_tc1_clean<-putative_tc1[,1:3]
write.table(putative_tc1_clean,file="putative_tc1_bayestyper.bed",row.names = FALSE,col.names = FALSE, quote=FALSE, sep="\t" )

### then it's possible to use bedtools to overlap these sequence with the genome using "getfasta"


## ----  upset plot  ----


library(UpSetR)

FE_SV<-read.table(file="Data/bayestyper_SV_info_table_Farmed_European.txt",header=TRUE)

WE_SV<-read.table(file="Data/bayestyper_SV_info_table_Wild_European.txt",header=TRUE)

WA_SV<-read.table(file="Data/bayestyper_SV_info_table_Wild_American.txt",header=TRUE)

FA_SV<-read.table(file="Data/bayestyper_SV_info_table_Farmed_American.txt",header=TRUE)




presence_FE<-ifelse(FE_SV$FREQ>0,1,0)
presence_FA<-ifelse(FA_SV$FREQ>0,1,0)
presence_WE<-ifelse(WE_SV$FREQ>0,1,0)
presence_WA<-ifelse(WA_SV$FREQ>0,1,0)

presence_FE<-ifelse(FE_SV$FREQ>0 & FE_SV$FREQ!=1,1,0)
presence_FA<-ifelse(FA_SV$FREQ>0 & FA_SV$FREQ!=1,1,0)
presence_WE<-ifelse(WE_SV$FREQ>0 & WE_SV$FREQ!=1,1,0)
presence_WA<-ifelse(WA_SV$FREQ>0 & WA_SV$FREQ!=1,1,0)


tableupset<-data.frame(CHROM=FE_SV$CHROM,POS=FE_SV$POS,END=FE_SV$END,Pres_FE=presence_FE,Pres_FA=presence_FA,Pres_WE=presence_WE,Pres_WA=presence_WA)

# Remove rows with NA values
tableupset <- na.omit(tableupset)


## make a table to see the counts

library(dplyr)

dfcount <- tableupset %>% group_by(Pres_FE,Pres_FA,Pres_WE,Pres_WA) %>% 
  summarise(total_count=n(),.groups = 'drop') %>%
  as.data.frame()

## make the plot

input <- c(
  Farmed_American = 8,
  Farmed_European = 67,
  Wild_European = 384,
  Wild_American = 16,
  "Wild_European&Farmed_European" = 1042,
  "Wild_European&Farmed_American" = 34,
  "Wild_European&Wild_American" = 65,
  "Farmed_American&Farmed_European" = 14,
  "Farmed_American&Wild_American" = 35,
  "Wild_American&Farmed_European" = 3,
  "Farmed_American&Wild_European&Wild_American" = 115,
  "Farmed_European&Wild_European&Wild_American" = 765,
  "Farmed_American&Farmed_European&Wild_American" = 10,
  "Farmed_American&Wild_European&Farmed_European" = 429,
  "Farmed_American&Wild_European&Wild_American&Farmed_European" = 2870
)



# Plot
upset(fromExpression(input), 
      nintersects = 40, 
      nsets = 6, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2.5, 
      point.size = 2.8, 
      line.size = 1
)



## ---- Sv overlapping annotations  ----

##import annotation
gffile<-read.table(file="Data/Salmo_salar.Ssal_v3.1.106_filtered.gff.gz",fill=TRUE)
### format the table
gffile<-gffile[,c(1,3,4,5,7)]
gffile<-gffile[gffile$V3%in%c("gene","mRNA","exon","CDS","pseudogene","lnc_RNA","snRNA","tRNA"),]
gffile$length<-as.numeric(as.character(gffile$V5))-as.numeric(as.character(gffile$V4))
gffile<-gffile[!is.na(gffile$length),]
gffile$ID<-paste(gffile$V1,gffile$V4,gffile$V5)



gffsummary <- gffile %>%
  group_by(V3) %>%  
  summarise(total_length = n())


### using bedtools, make an intersect between the annotation and the SVs coordinate

## then import the results

annotfile<-read.table(file="Data/overlap_SV_annot_v2.txt",fill=TRUE,sep="\t",quote = "")
## and format the table
annotfile<-annotfile %>% distinct(V1, V3,V4,V5, .keep_all = TRUE)
annotfile$IDSV<-paste(annotfile$V10,annotfile$V11,annotfile$V12,sep="_")
annotfile<-annotfile[,c(1,3,4,5,7,13)]
annotfile<-annotfile[annotfile$V3%in%c("gene","mRNA","exon","CDS","pseudogene","lnc_RNA","snRNA","tRNA"),]


annotfile$V4<-as.numeric(as.character(annotfile$V4))
annotfile$V5<-as.numeric(as.character(annotfile$V5))
annotfile$length<-annotfile$V5-annotfile$V4
annotfile<-annotfile[!is.na(annotfile$length),]
annotfile$ID<-paste(annotfile$V1,annotfile$V4,annotfile$V5)


overlapsummary <- annotfile %>%
  group_by(V3) %>%  
  summarise(total_length = n())



summarytotal<-data.frame(type=gffsummary$V3,nboverlap=overlapsummary$total_length,nbannot=gffsummary$total_length)
### make a ratio of the the number of annotations overlapped by Sv / the number of total annotation in the genome
summarytotal$ratiooverlap<-summarytotal$nboverlap/summarytotal$nbannot


### plot results
ggplot(summarytotal,aes(x=type,y=ratiooverlap,fill=type))+geom_bar(stat="identity")+theme_bw(15)+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+geom_text(aes(label = signif(nboverlap)), nudge_y = 0.003, size=5)+
  scale_fill_brewer(palette = "Set2")