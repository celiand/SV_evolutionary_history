##### Celian Diblasi
#### Transposable elements overlap with SV
#### 28/10/2025

#### to run the distribution for random Svs, use repeatSvanalysis.sh with this R script

## ---- First we need to make the two sets of data : real SVs and random SVs, with their location ----

### For real Sv , just extract SV boundaries
FE_SV_BT<-read.table(file="Data/bayestyper_SV_info_table_Farmed_European.txt",header=TRUE)
FE_SV_BT$ID<-paste(FE_SV_BT$CHROM,FE_SV_BT$POS,FE_SV_BT$END,sep="_")
FE_SV_BT$length<-FE_SV_BT$END-FE_SV_BT$POS
FE_SV_BT<-FE_SV_BT[FE_SV_BT$length>=50,]

SV_boundaries_table<-data.frame(CHROM=c(FE_SV_BT$CHROM),START=c(FE_SV_BT$POS),END=c(FE_SV_BT$END),IDSV=c(FE_SV_BT$ID))

options("scipen"=100, "digits"=4) ##disable scientific notation

write.table(SV_boundaries_table,file="SV_boundaries_table.bed",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")



## For random Svs, just replace each Sv randomly
### (call this script in a loop from an external script to make different random SV files)
FE_SV_BT<-read.table(file="Data/bayestyper_SV_info_table_Farmed_European.txt",header=TRUE)
FE_SV_BT$ID<-paste(FE_SV_BT$CHROM,FE_SV_BT$POS,FE_SV_BT$END,sep="_")
FE_SV_BT$length<-FE_SV_BT$END-FE_SV_BT$POS
FE_SV_BT<-FE_SV_BT[FE_SV_BT$length>=50,]

## get the chromosome names
chromosomes<-unique(FE_SV_BT$CHROM)
### and the chromosomes size
lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)

FE_SV_BT$chromnumer<-as.numeric(as.character(ifelse(substr(FE_SV_BT$CHROM,4,4)==0,substr(FE_SV_BT$CHROM,5,5),substr(FE_SV_BT$CHROM,4,5))))


for (i in 1:length(FE_SV_BT$CHROM)){
  lengthSV<-FE_SV_BT[i,"length"]
  SVplaced<-FALSE
  randomchrom<-FE_SV_BT[i,"chromnumer"]
  while(SVplaced==FALSE){
    ### randomly place the Sv in the chromosome
    randompos<-sample(10000:lengthvector[randomchrom],1)
    randomend<-randompos+lengthSV
    if(randomend<=lengthvector[randomchrom]-10000){
      ### just check if the start of the Sv is not too close to the end of the chromosome
      SVplaced<-TRUE
    }
  }
  FE_SV_BT[i,"falseSVchrom"]<-chromosomes[randomchrom]
  FE_SV_BT[i,"falseSVpos"]<-randompos
  FE_SV_BT[i,"falseSVend"]<-randompos+lengthSV
}

RandomSVs<-data.frame(CHROM=FE_SV_BT$falseSVchrom,start=FE_SV_BT$falseSVpos,end=FE_SV_BT$falseSVend)
### save the table
write.table(RandomSVs,file="Data/RandomSVs_bis.bed",row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")

randomssv<-read.table(file="Data/RandomSVs_bis.bed")
randomssv$ID<-paste(randomssv$V1,randomssv$V2,randomssv$V3,sep="_")

## extract boudnaries as for real observed SVs and register them
random_SV_boundaries_table<-data.frame(CHROM=c(randomssv$V1),START=c(randomssv$V2),END=c(randomssv$V3),IDSV=c(randomssv$ID))

### use this file to do the following steps
write.table(random_SV_boundaries_table,file="Data/random_SV_boundaries_table.bed",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")

##### note: for flaking regions of SVs, use the part below
extended_boundaries<-1000

random_SV_extended_boundaries_before<-data.frame(CHROM=c(randomssv$V1),START=c(randomssv$V2-extended_boundaries),END=c(randomssv$V2),IDSV=c(randomssv$ID),Location="before")
random_SV_extended_boundaries_after<-data.frame(CHROM=c(randomssv$V1),START=c(randomssv$V3),END=c(randomssv$V3+extended_boundaries),IDSV=c(randomssv$ID),Location="after")

random_SV_extended_boundaries_table<-rbind.data.frame(random_SV_extended_boundaries_before,random_SV_extended_boundaries_after)

write.table(random_SV_extended_boundaries_table,file="Data/random_SV_extended_boundaries_table.bed",row.names = FALSE,quote=FALSE,col.names = FALSE,sep="\t")





## ---- then we need to initialize the tables with the real observed data ---- 
### import overlap of observed Sv with TEs
real_overlap<-read.table(file="Data/SV_boundaries_overlap_repeat_realdata.txt", sep="\t")

nb_sv_overlap<-length(unique(real_overlap$V4))

nb_sv_overlap_SDonly<-length(unique(real_overlap[real_overlap$V8=="Segmental duplication","V4"]))

nb_sv_overlap_nosd<-length(unique(real_overlap[real_overlap$V8!="Segmental duplication","V4"]))

nb_sv_overlap_LINE<-length(unique(real_overlap[grepl("^LINE", real_overlap$V8),"V4"]))

nb_sv_overlap_DNA<-length(unique(real_overlap[grepl("^DNA", real_overlap$V8),"V4"]))

nb_sv_overlap_SINE<-length(unique(real_overlap[grepl("^SINE", real_overlap$V8),"V4"]))

nb_sv_overlap_LTR<-length(unique(real_overlap[grepl("^LTR", real_overlap$V8),"V4"]))

new_df_nb<-data.frame(condition="real",typerepeat=c("All","SD","NoSD","LINE","SINE","DNA","LTR"),nb_of_sv=c(nb_sv_overlap,nb_sv_overlap_SDonly,nb_sv_overlap_nosd,nb_sv_overlap_LINE,nb_sv_overlap_SINE,nb_sv_overlap_DNA,nb_sv_overlap_LTR),runocc=0)

### this table will be used to add data from random SVs
write.table(new_df_nb,file="df_nb_sv_repeat_global.txt",row.names = FALSE)



#### ---- Use this to compute proportion of overlap (within a loop) ----

args <- commandArgs()


run<- args[6]


random_overlap<-read.table(file="Data/random_SV_boundaries_table_overlap_repeat.txt", sep="\t")

nb_sv_overlap<-length(unique(random_overlap$V4))

nb_sv_overlap_SDonly<-length(unique(random_overlap[random_overlap$V8=="Segmental duplication","V4"]))

nb_sv_overlap_nosd<-length(unique(random_overlap[random_overlap$V8!="Segmental duplication","V4"]))

nb_sv_overlap_LINE<-length(unique(random_overlap[grepl("^LINE", random_overlap$V8),"V4"]))

nb_sv_overlap_DNA<-length(unique(random_overlap[grepl("^DNA", random_overlap$V8),"V4"]))

nb_sv_overlap_SINE<-length(unique(random_overlap[grepl("^SINE", random_overlap$V8),"V4"]))

nb_sv_overlap_LTR<-length(unique(random_overlap[grepl("^LTR", random_overlap$V8),"V4"]))

new_df_nb<-data.frame(condition="random",typerepeat=c("All","SD","NoSD","LINE","SINE","DNA","LTR"),nb_of_sv=c(nb_sv_overlap,nb_sv_overlap_SDonly,nb_sv_overlap_nosd,nb_sv_overlap_LINE,nb_sv_overlap_SINE,nb_sv_overlap_DNA,nb_sv_overlap_LTR),runocc=run)

prev_table_stat<-read.table(file="Data/df_nb_sv_repeat_global.txt",header=TRUE)

new_table_stat_complete<-rbind.data.frame(prev_table_stat,new_df_nb,stringsAsFactors = FALSE)

write.table(new_table_stat_complete,file="Data/df_nb_sv_repeat_global.txt",row.names = FALSE)


#### use the same logic with the flaking region to make "Data/df_nb_sv_repeat_global_extended.txt"

#### ---- plot repeat overlap  ----

#### first within SV

table_stat_full<-read.table(file="Data/df_nb_sv_repeat_global.txt",header=TRUE)
table_stat_full$proportion<-table_stat_full$nb_of_sv/4244

library(ggplot2)

factororder <- c(table_stat_full[1:7,2])
table_stat_full$typerepeat <- factor(table_stat_full$typerepeat , levels = factororder)

realdata_sub<-table_stat_full[table_stat_full$condition=="real" & table_stat_full$typerepeat!="NoSD",]

plotdata<-table_stat_full[table_stat_full$condition=="random" & table_stat_full$typerepeat!="NoSD",]


ggplot(plotdata, aes(x = factor(typerepeat, levels = rev(unique(typerepeat))), y = proportion)) +
  geom_violin(alpha = 0.6,fill="lightblue",scale = "width") + 
  geom_point(data = realdata_sub, aes(x = factor(typerepeat, levels = rev(unique(typerepeat))), y = proportion), 
             color = "firebrick", size = 3, shape = 18) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(0,1)+theme_bw(15)


### then in flanking regions


table_stat_full_out<-read.table(file="Data/df_nb_sv_repeat_global_extended.txt",header=TRUE)
table_stat_full_out$proportion<-table_stat_full_out$nb_of_sv/4244


factororder <- c(table_stat_full_out[1:7,2])
table_stat_full_out$typerepeat <- factor(table_stat_full_out$typerepeat , levels = factororder)

realdata_sub_out<-table_stat_full_out[table_stat_full_out$condition=="real" & table_stat_full_out$typerepeat!="NoSD",]

plotdata_out<-table_stat_full_out[table_stat_full_out$condition=="random" & table_stat_full_out$typerepeat!="NoSD",]




ggplot(plotdata_out, aes(x = factor(typerepeat, levels = rev(unique(typerepeat))), y = proportion)) +
  geom_violin(alpha = 0.6,fill="lightblue",scale = "width") + 
  geom_point(data = realdata_sub_out, aes(x = factor(typerepeat, levels = rev(unique(typerepeat))), y = proportion), 
             color = "firebrick", size = 4, shape = 18) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylim(0,1)+theme_bw(15)




#### ---- script used to make a common files for all repeats ----

### chromosomes length
lengthvector<-c(174498729,95481959,105780080,90536438,92788608,96060288,68862998,28860523,161282225,125877811,111868677,101677876,
                114417674,101980477,110670232,96486271,87489397,84084598,88107222,96847506,59819933,63823863,52460201,49354470,54385492,
                55994222,45305548,41468476,43051128)

## TE library (see methods)
maskeddata<-read.table(file="/mnt/SCRATCH/cedi/phDSalmon/TE_library/masked_genome/Salmo_salar.Ssal_v3.1.dna_sm.toplevel.fa.out",skip=1,header=TRUE,fill=TRUE,row.names = NULL)
## format tables
maskeddata<-maskeddata[,c(2,5,6,7,11)]
colnames(maskeddata)<-c("DivFromCons","chromgen","startgen","endgen","family")

### get segmental duplication
##from https://genome-euro.ucsc.edu/cgi-bin/hgTables?hgsid=342142508_Q89ACLALF2Efu55Yz3VPYRaXf[b