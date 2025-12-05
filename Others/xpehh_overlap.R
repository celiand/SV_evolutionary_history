### check overlap of SVs with xpehh data from buso et al.



# ---- make separate analysis for EU and AM ---- 
## ---- EU ----

### read global data

allxpehh<-read.table(file="Data/xpehh_table_EU_3_pau_formatted.txt")

## format SNp position so it's on a continuous scale
lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector+50000000)[1:28])

chromname<-unique(allxpehh$V1)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
allxpehh$newposSNP<-allxpehh$V2+listval[match(allxpehh$V1,listval$chromname),3] 

## same for lostruct SVs overlapped with xpehh (made using bedtools)
euoverlap<-read.table(file="Data/PCA_SV_xpehh_EU.txt")

euoverlap$newposSNP<-euoverlap$V5+listval[match(euoverlap$V1,listval$chromname),3] 

## same for bayestyper SVs
euoverlapBT<-read.table(file="Data/BT_SV_xpehh_EU.txt")

euoverlapBT$newposSNP<-euoverlapBT$V5+listval[match(euoverlapBT$V1,listval$chromname),3] 

allxpehh$type<-ifelse(allxpehh$newposSNP%in%euoverlap$newposSNP,ifelse(allxpehh$newposSNP%in%euoverlapBT$newposSNP,"both","PCA_SV"),ifelse(allxpehh$newposSNP%in%euoverlapBT$newposSNP,"BT_SV","background"))


library(ggplot2)

##plot 

ggplot(allxpehh, aes(x = newposSNP, y = V5)) +
  geom_point(data = subset(allxpehh, type == "background"), aes(color = type, size = type, alpha = type)) +
  geom_point(data = subset(allxpehh, type == "PCA_SV"), aes(color = type, size = type, alpha = type)) +
  geom_point(data = subset(allxpehh, type == "BT_SV"), aes(color = type, size = type, alpha = type)) +
  scale_size_manual(values = c(0.2, 1.5, 1.5)) +
  scale_alpha_manual(values = c(0.2, 1,1)) +
  theme_minimal()+scale_color_manual(values=c("grey70","dodgerblue3","firebrick"))+ylab("Log pvalue")+
  ggtitle("European XPEHH overlap")+theme(axis.text.x = element_text(angle = 45, hjust=1))+
  scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+geom_hline(yintercept = -log10(0.0001),color="black",linewidth=1, linetype="dashed")

allxpehhsub<-allxpehh[allxpehh$V5>4 & allxpehh$type!="background",]


### make a table with region with significant signal
data_SV_region_significant<-data.frame(chrom=c(allxpehhsub[c(1,2,3),"V1"]),start=c(allxpehhsub[c(1,2,3),"V2"]), end= c(allxpehhsub[c(1,2,69),"V3"]), typeSV=c(allxpehhsub[c(1,2,3),"type"]))
data_SV_region_significant$ID<-paste(data_SV_region_significant$chrom,data_SV_region_significant$start,data_SV_region_significant$end,sep="_")
data_SV_region_significant$region<-"EU"




## same analysis for AM 
## ---- AM ----
### read global data

allxpehh<-read.table(file="Data/xpehh_table_AM_3_pau_formatted.txt")

lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)
chromadd<-c(0,cumsum(lengthvector+50000000)[1:28])

chromname<-unique(allxpehh$V1)
chrommid<-chromadd+lengthvector/2
listval<-cbind.data.frame(chromname,lengthvector,chromadd,chrommid)
allxpehh$newposSNP<-allxpehh$V2+listval[match(allxpehh$V1,listval$chromname),3] 


euoverlap<-read.table(file="Data/PCA_SV_xpehh_AM.txt")

euoverlap$newposSNP<-euoverlap$V5+listval[match(euoverlap$V1,listval$chromname),3] 


euoverlapBT<-read.table(file="Data/BT_SV_xpehh_AM.txt")

euoverlapBT$newposSNP<-euoverlapBT$V5+listval[match(euoverlapBT$V1,listval$chromname),3] 

allxpehh$type<-ifelse(allxpehh$newposSNP%in%euoverlap$newposSNP,ifelse(allxpehh$newposSNP%in%euoverlapBT$newposSNP,"both","PCA_SV"),ifelse(allxpehh$newposSNP%in%euoverlapBT$newposSNP,"BT_SV","background"))


library(ggplot2)


ggplot(allxpehh, aes(x = newposSNP, y = V5)) +
  geom_point(data = subset(allxpehh, type == "background"), aes(color = type, size = type, alpha = type)) +
  geom_point(data = subset(allxpehh, type == "PCA_SV"), aes(color = type, size = type, alpha = type)) +
  geom_point(data = subset(allxpehh, type == "BT_SV"), aes(color = type, size = type, alpha = type)) +
  scale_size_manual(values = c(0.2, 1.5, 1.5)) +
  scale_alpha_manual(values = c(0.2, 1,1)) +
  theme_minimal()+scale_color_manual(values=c("grey70","dodgerblue3","firebrick"))+
  ggtitle("American XPEHH overlap")+theme(axis.text.x = element_text(angle = 45, hjust=1))+ylab("Log pvalue")+
  scale_x_continuous(label =listval$chromname, breaks = listval$chrommid )+geom_hline(yintercept = -log10(0.0001),color="black",linewidth=1, linetype="dashed")


allxpehhsub<-allxpehh[allxpehh$V5>4 & allxpehh$type!="background",]


## also extract boundaries of region with positive signals
data_SV_region_significant_AM<-data.frame(chrom=c("ssa04","ssa13","ssa16"),start=c(29044783,76644638,60662249), end= c(29238337,76691235,60879497), typeSV=c("PCA_SV","PCA_SV","PCA_SV"))
data_SV_region_significant_AM$ID<-paste(data_SV_region_significant_AM$chrom,data_SV_region_significant_AM$start,data_SV_region_significant_AM$end,sep="_")
data_SV_region_significant_AM$region<-"AM"


# ----  save significant regions ----

data_SV_region_significant<-rbind.data.frame(data_SV_region_significant,data_SV_region_significant_AM,stringsAsFactors = FALSE)

##change chromosome name to be compatible with gff
data_SV_region_significant <- data_SV_region_significant %>%
  mutate(chrom = as.numeric(gsub("ssa", "", chrom)))

options("scipen"=100, "digits"=4) ##disable scientific notation

write.table(data_SV_region_significant,file="SV_overlap_xpehh.bed",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")


### check which genes are overlapped (after overlap with bedtools and the genome annotation)
overlaptable<-read.table(file="Data/SV_overlap_xpehh_genes.txt",header=FALSE,sep="\t")
overlaptable_gene<-overlaptable[overlaptable$V9=="gene",]

write.table(overlaptable_gene,file="gene_overlap_xpehh.txt",row.names = FALSE,col.names = FALSE,sep="\t")