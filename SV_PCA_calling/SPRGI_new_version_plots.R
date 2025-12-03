## first get data from the script
options("scipen"=100, "digits"=4) ##disable scientific notation

args <- commandArgs()


region <- args[6]

##region<-"ssa05_44020000_44084228" ###if you want to run a single region

chromosome<-sapply(strsplit(region,"_"),FUN = `[[`, 1)
start<-as.numeric(as.character(sapply(strsplit(region,"_"),FUN = `[[`, 2)))
end<-as.numeric(as.character(sapply(strsplit(region,"_"),FUN = `[[`, 3)))



regionfile<-read.table(file="data/interesting_region_bis.txt",header=TRUE)


coderegion<-regionfile[regionfile$CHROM==chromosome & regionfile$Manual_start==start,"newcode"]


diff<-end-start

sizetochange<-0.5*diff

newstart<-start-sizetochange
newend<-end+sizetochange

newregion<-paste(chromosome,newstart,newend,sep="_")




###plot genotype
region_data_raw_filename<-paste("data/sprgi_vcffiles/PCA_",newregion,"_int.vcf.gz_258.stat",sep="")

##in case the file does not exists, get the correct name
if(!file.exists(region_data_raw_filename)){
  newregion<-paste(chromosome,floor(newstart),floor(newend),sep="_")
  region_data_raw_filename<-paste("data/sprgi_vcffiles/PCA_",newregion,"_int.vcf.gz_258.stat",sep="")
  
}

if(!file.exists(region_data_raw_filename)){
  newregion<-paste(chromosome,floor(newstart),ceiling(newend),sep="_")
  region_data_raw_filename<-paste("data/sprgi_vcffiles/PCA_",newregion,"_int.vcf.gz_258.stat",sep="")
  
}

if(!file.exists(region_data_raw_filename)){
  newregion<-paste(chromosome,ceiling(newstart),floor(newend),sep="_")
  region_data_raw_filename<-paste("data/sprgi_vcffiles/PCA_",newregion,"_int.vcf.gz_258.stat",sep="")
  
}

if(!file.exists(region_data_raw_filename)){
  newregion<-paste(chromosome,ceiling(newstart),ceiling(newend),sep="_")
  region_data_raw_filename<-paste("data/sprgi_vcffiles/PCA_",newregion,"_int.vcf.gz_258.stat",sep="")
  
}

region_data_raw<-read.table(file=region_data_raw_filename,header=FALSE)

filenamepcagt<-paste("data/sprgi_gt_pca_files/",region,".txt",sep="")

region_genotype<-read.table(file=filenamepcagt,header=TRUE)

library("gplots")
library(ggplot2); library(reshape2)
library(tidyverse)

vcf1<-region_data_raw

vcf1[vcf1=="0/0"]<- "0"
vcf1[vcf1=="0/1"]<- "1"
vcf1[vcf1=="1/0"]<- "1"
vcf1[vcf1=="1/1"]<- "2"
vcf1[vcf1=="./."]<- "0"
vcf1[vcf1=="1|0"]<- "1"
vcf1[vcf1=="0|1"]<- "1"
vcf1[vcf1=="1|1"]<- "2"


vcf3 <- data.frame(apply(vcf1, 2, function(x) as.numeric(as.character(x))))

vcf2 <- t(vcf3)
vcf4 <- data.matrix(vcf2)
colnames(vcf4) <- factor(row.names(vcf1))

vcf4 %>% janitor::row_to_names(1) -> vcf4

##### haplotype plot

my_palette <- colorRampPalette(c("#6b562e","#a86f03","#d8c194")) (n=3)

region_genotype$colorsGT <- ifelse(region_genotype$GT == "NI/NI", "#baba09",
                                   ifelse(region_genotype$GT== "NI/I", "#d22121",
                                          ifelse(region_genotype$GT == "I/I", "#3274e4", NA)))


plot_heatmap_filename2<-paste("data/sprgi_plots/",coderegion,"_heatmap",".png",sep="") ##filename under which to register the plot

colorspop<-region_genotype$colorsGT

png(plot_heatmap_filename2, width=900, height=600)

colnamesvcf<-colnames(vcf4)
numericcoord<-as.numeric(as.character(colnamesvcf))
###define sv and flanking regions
flanking_a<-length(numericcoord[numericcoord<start])
flanking_b<-length(numericcoord[numericcoord>end])
SV_region<-length(numericcoord[numericcoord>=start & numericcoord<=end])
col_SV_leg<-c(rep("grey70",flanking_a),rep("black",SV_region),rep("grey70",flanking_b))

###plot
heatmap.2(vcf4,trace="none", na.color = "#e7e7e7",scale="none", margins=c(12,8),
          col = my_palette,density.info="none",dendrogram = "none", RowSideColors=colorspop,Rowv=TRUE, Colv=FALSE,main=region,ColSideColors=col_SV_leg)
dev.off()



### same haplotype plot but for population color


region_genotype$colorspop <- ifelse(region_genotype$pop == "Farmed_European", "#1f78b4",
                                   ifelse(region_genotype$pop== "Farmed_American", "#33a02c",
                                          ifelse(region_genotype$pop == "Wild_European", "#a6cee3",  
                                                 ifelse(region_genotype$pop == "Wild_American", "#b2df8a", NA))))

plot_heatmap_filename3<-paste("data/sprgi_plots/",coderegion,"_heatmap_pop",".png",sep="")

colorspop<-region_genotype$colorspop

png(plot_heatmap_filename3, width=900, height=600)
heatmap.2(vcf4, trace="none", na.color = "#e7e7e7",scale="none", margins=c(12,8),
          col = my_palette,density.info="none",dendrogram = "none", RowSideColors=colorspop,Rowv=TRUE, Colv=FALSE,main=region,ColSideColors=col_SV_leg)
dev.off()




### plot heterozygosity
library("adegenet")
library("hierfstat")
library("pegas")
library(vcfR)
library(ggplot2)

region_data_raw<-read.table(file=region_data_raw_filename,header=FALSE)

filename_vcf<-paste("data/sprgi_vcffiles/PCA_",newregion,"_int.vcf.gz",sep="")

vcf_data<-read.vcfR(file=filename_vcf)

genind_obj<-vcfR2genind(vcf_data)

pop(genind_obj)<-region_genotype$GT

## data for each genotype
div_II <- summary(genind_obj[genind_obj$pop=="I/I"])

data_II <- data.frame(
  Loci_Position = region_data_raw$V1,  
  Observed_Heterozygosity = div_II$Hobs,
  Genotype="I/I"
)



div_NII <- summary(genind_obj[genind_obj$pop=="NI/I"])

data_NII <- data.frame(
  Loci_Position = region_data_raw$V1,  
  Observed_Heterozygosity = div_NII$Hobs,
  Genotype="NI/I"
)

div_NINI <- summary(genind_obj[genind_obj$pop=="NI/NI"])

data_NINI <- data.frame(
  Loci_Position = region_data_raw$V1,  
  Observed_Heterozygosity = div_NINI$Hobs,
  Genotype="NI/NI"
)

###merge all genotype data
data_het_full<-rbind.data.frame(data_II,data_NII,data_NINI,stringsAsFactors = FALSE)

###general plot
plothet<-ggplot(data_het_full, aes(x = Loci_Position, y = Observed_Heterozygosity,color=Genotype, shape=Genotype)) +
  geom_point(alpha=0.5,size=3.5) +
  labs(x = "Loci Position", y = "Observed Heterozygosity")+scale_color_manual(values=c("#3274e4","#d22121","#baba09"))+theme_bw(30)+xlim(c(newstart,newend))+annotate("rect",xmin = newstart, xmax = start,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45) +annotate("rect",xmin = end, xmax = newend,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45)#+ geom_vline(xintercept = c(start, end), color = "black",linetype="dashed",size=2.5,alpha=0.5)




###plot sprgi and dxy
pidatafile<-paste("data/sprgi_pidxy/pixy_pi_",newregion,".txt",sep="")
pidata<-read.table(file=pidatafile,header=TRUE)
pidatafilter<-pidata[pidata$no_sites>=10,]

piplot<-ggplot(pidatafilter,aes(x=window_pos_1,y=avg_pi,color=pop))+scale_linetype_manual(values = c("solid", "dashed", "dotted"))+
  scale_shape_manual(values = c(16, 17, 15))+scale_color_manual(values=c("#3274e4","#d22121","#baba09"))+theme_bw(30)+theme(axis.text.x = element_text(angle = 45, hjust=1))+ annotate("rect",xmin = newstart, xmax = start,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45) +annotate("rect",xmin = end, xmax = newend,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45)+xlim(c(newstart,newend))+geom_line(aes(linetype = pop, color = pop),size=2.5)#+geom_vline(xintercept = c(start, end), color = "blue",linetype="dashed",size=2.5)


dxydatafile<-paste("data/sprgi_pidxy/pixy_dxy_",newregion,".txt",sep="")
dxydata<-read.table(file=dxydatafile,header=TRUE)
dxydatafilter<-dxydata[dxydata$no_sites>=10,]

dxydatafilter$comparison<-paste(dxydatafilter$pop1,dxydatafilter$pop2,sep=" vs ")

dxyplot<-ggplot(dxydatafilter,aes(x=window_pos_1,y=avg_dxy,color=comparison))+scale_linetype_manual(values = c("solid", "dashed", "dotted"))+
  +theme_bw(30)+theme(axis.text.x = element_text(angle = 45, hjust=1))+scale_color_manual(values=c("#c064be","#539034","#ec9f23"))+xlim(c(newstart,newend))+geom_line(aes(linetype = comparison, color = comparison),size=2.5)+ annotate("rect",xmin = newstart, xmax = start,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45) +annotate("rect",xmin = end, xmax = newend,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45)#+ geom_vline(xintercept = c(start, end), color = "blue",linetype="dashed",size=2.5)


###plot fst between each genotype

fst_filename_1<-paste("data/sprgi_fstfiles/",newregion,"_fst_II_vs_NINI.txt.weir.fst",sep="")
fst_region_1<-read.table(file=fst_filename_1,header=TRUE)
fst_region_1$comparison="II_vs_NINI"

fst_filename_2<-paste("data/sprgi_fstfiles/",newregion,"_fst_II_vs_NII.txt.weir.fst",sep="")
fst_region_2<-read.table(file=fst_filename_2,header=TRUE)
fst_region_2$comparison="II_vs_NII"

fst_filename_3<-paste("data/sprgi_fstfiles/",newregion,"_fst_NII_vs_NINI.txt.weir.fst",sep="")
fst_region_3<-read.table(file=fst_filename_3,header=TRUE)
fst_region_3$comparison="NII_vs_NINI"

##merge all fst data
fst_full_data<-rbind.data.frame(fst_region_1,fst_region_2,fst_region_3,stringsAsFactors = FALSE)

library(ggplot2)
##plot
plot_FST<-ggplot(fst_full_data,aes(x=POS,y=WEIR_AND_COCKERHAM_FST,color=comparison,shape=comparison))+
  scale_shape_manual(values = c(16, 17, 15))+geom_point(alpha=0.6,size=3.5)+theme_bw(30)+scale_color_manual(values=c("#c064be","#539034","#ec9f23"))+xlim(c(newstart,newend))+ annotate("rect",xmin = newstart, xmax = start,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45) +annotate("rect",xmin = end, xmax = newend,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45)#+ geom_vline(xintercept = c(start, end), color = "blue",linetype="dashed",size=2.5)



### add depth plot for the SV

library(tidyverse)
library(ggplot2)
depthdata<-read.table(file="data/subdepthdatassa03.txt",header=TRUE) ## extract data for the region
regionfile<-read.table(file="data/interesting_region_bis.txt",header=TRUE)

i<-3 ## number of the region


region_genotype<-read.table(file=filenamepcagt,header=TRUE)
depthpart<-depthdata[depthdata$chrom==chromosome & depthdata$start>=newstart & depthdata$end<=newend, ]
depthpart$sample <- gsub("_RG", "", depthpart$sample)
merged_depthpart <- depthpart %>%
  left_join(region_genotype, by = c("sample" = "IID"))
pdepth<-ggplot(merged_depthpart, aes(x = start, y = rrd, color = GT, group = interaction(start, GT))) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal(base_size = 11)+
  labs(
    x = "Genomic window start", y = "RRD"
  ) +
  scale_fill_brewer(palette = "Set2")+theme_bw(20)+scale_color_manual(values=c("#3274e4","#d22121","#baba09"))+
  annotate("rect",xmin = newstart, xmax = start,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45) +annotate("rect",xmin = end, xmax = newend,ymin = -Inf, ymax = Inf, fill = "grey", alpha = 0.45)+xlim(c(newstart,newend))



### combine the plots together (heterozygosity, pi, dxy, FST, depth)
library(cowplot)



# Combine the plots
combined_plot <- plot_grid(plot_FST,dxyplot,piplot,plothet,pdepth, ncol = 1, align = "v")
combinedplotname<-paste(coderegion,"_completeplot",".png",sep="")

###save the plot
ggsave(filename=combinedplotname,plot=combined_plot,device="png",path="/mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/pauline_data/sprgi_plots/",width=30,height=34,units="in",limitsize = FALSE)






## then make the pca plot

library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)



## retrieve full name of file


region_pca_filename_eigenvec<-paste("data/sprgi_pcafiles/PCA_",region,".eigenvec",sep="")
region_pca_filename_eigenval<-paste("data/sprgi_pcafiles/PCA_",region,".eigenval",sep="")

##then plot
region_pca <-data.table::fread(region_pca_filename_eigenvec[1])  #for entire PCA plot

region_pca_var <-data.table::fread(region_pca_filename_eigenval[1]) 
colnames(region_pca_var)<-"eigenval"
region_pca_var<- region_pca_var %>%
  mutate(var.prop= round(eigenval/sum(eigenval),2),PC=1:n()) 

prop2<-as.numeric(as.character(region_pca_var[2,2]))*100
prop1<-as.numeric(as.character(region_pca_var[1,2]))*100

#input the populations file
popfilename<-paste("data/sprgi_gt_pca_files/",region,".txt",sep="")
pops_file <- read.table(popfilename,header=TRUE)

#merge the two tables to have location info on the pca + legend by pop or type of samples
PCA_file_merged <- merge(region_pca, pops_file, by.x = "#IID",by.y = "IID")

#Plot

ylabel<-paste("PC2 - ",prop2,"%",sep="")
xlabel<-paste("PC1 - ",prop1,"%",sep="")

plotPCA<-ggplot(PCA_file_merged,aes(x = PC1, y = PC2, group=GT)) + geom_point((aes(colour = GT, shape=pop)), size = 4, show.legend = TRUE) + 
  theme_bw(25) + labs(x = xlabel, y = ylabel) +
  scale_color_manual(values = c("#3274e4","#d22121","#baba09"))+scale_shape_manual(values=c(16, 15, 1,0))


plot_pca_filename<-paste(coderegion,"_PCA",region,".png",sep="")

ggsave(filename=plot_pca_filename,plot=plotPCA,device="png",path="/mnt/SCRATCH/cedi/phDSalmon/evolutionary_analysis/pauline_data/sprgi_plots/",width=30,height=24,units="in")



