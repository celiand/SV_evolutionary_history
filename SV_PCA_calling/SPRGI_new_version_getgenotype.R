## first get data from the script
options("scipen"=100, "digits"=4) ##disable scientific notation

args <- commandArgs()


region <- args[6]

#region<-"ssa01_92700000_92900000"

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

###let's get the genotype

library(invClust)
library(snpStats)



region_data_raw_filename<-paste("data/sprgi_vcffiles/PCA_",region,"_int.vcf.gz_258.stat",sep="")
region_data_raw<-read.table(file=region_data_raw_filename,header=FALSE)
samplefile<-read.table(file="data/samples_int_region_bis.txt",header=FALSE)
sampname<-samplefile$V1
colnames(region_data_raw)<-c("POS",sampname)

vcf1<-region_data_raw

vcf1[vcf1=="0/0"]<- 1
vcf1[vcf1=="0/1"]<- 2
vcf1[vcf1=="1/0"]<- 2
vcf1[vcf1=="1/1"]<- 3
vcf1[vcf1=="1|0"]<- 2
vcf1[vcf1=="0|1"]<- 2
vcf1[vcf1=="1|1"]<- 3
vcf1[vcf1=="0|0"]<- 1
vcf1[vcf1=="./."]<- 0



testmat<-t(vcf1)
colnames(testmat)<-testmat[1,]
testmat<-testmat[-c(1),]
region_snpmatrix<-new("SnpMatrix",testmat)

annotfile_name<-paste("data/sprgi_mapfiles/mapfile_",region,".map",sep="")
annotfile<-read.table(file=annotfile_name)

annotfile2<-data.frame(chromosome=annotfile$V1,snp.name=annotfile$V4,position=annotfile$V4)

### investigate region of interest
roi2<-data.frame(chr=chromosome,LBP=start,RBP=end,reg=region)

svint<-invClust(roi=roi2,wh=1,geno=region_snpmatrix,annot=annotfile2,dim=1)

genotype<-as.data.frame(svint["genotypes"])






get_max_column <- function(row) {
  col_names <- names(row)
  max_col_index <- which.max(row)
  return(col_names[max_col_index])
}

# Apply the function to each row of the data frame
max_columns <- apply(genotype, 1, get_max_column)

populations<-c(rep("Farmed_European", 112),rep("Farmed_American",80),rep("Wild_European",98),rep("Wild_American",79))

region_genotype<-data.frame(IID=c(row.names(genotype)),GT=max_columns,pop=populations)

filenamepcagt<-paste("data/sprgi_gt_pca_files/",region,".txt",sep="")
write.table(region_genotype,file=filenamepcagt,row.names = FALSE)

region_genotype_nopop<-data.frame(IID=c(row.names(genotype)),GT=max_columns)
filenamepcagtmopop<-paste("data/sprgi_gt_pca_files/gtnopop_",region,".txt",sep="")
write.table(region_genotype_nopop,file=filenamepcagtmopop,row.names = FALSE,sep="\t",quote=FALSE,col.names = FALSE)


##register individuals in separate files for their genotypes
II_filename<-paste("data/sprgi_gt_AABB/",region,"_II_ind.txt",sep="")
NINI_filename<-paste("data/sprgi_gt_AABB/",region,"_NINI_ind.txt",sep="")
NII_filename<-paste("data/sprgi_gt_AABB/",region,"_NII_ind.txt",sep="")

IIind<-region_genotype[region_genotype$GT=="I/I","IID"]
NINIind<-region_genotype[region_genotype$GT=="NI/NI","IID"]
NIIind<-region_genotype[region_genotype$GT=="NI/I","IID"]

write.table(IIind,file=II_filename,row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(NINIind,file=NINI_filename,row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(NIIind,file=NII_filename,row.names = FALSE,col.names = FALSE,quote = FALSE)


