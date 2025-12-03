##read arguments and get name of region processed
args <- commandArgs()


name <- sub("\"(.*)\"", "\\1", args[6])


alldata<-read.table(file="data/full_info_windows_PCA.txt",header=TRUE)

regionprocessed<-alldata[alldata$code==name,]



## extract individuals

for( i in 1:max(regionprocessed$cluster)){
  indiv<-regionprocessed[regionprocessed$cluster==i,"ind"]
  namefile1<-paste(name,i,sep="_")
  namefile2<-paste("data/tempfilefst/",namefile1,".txt",sep="")
  write.table(indiv,file=namefile2,quote=FALSE,row.names=FALSE,col.names = FALSE)
}



## extract the region coordinates

coordinates<-data.frame(chrom=regionprocessed[1,"chrom"],start=regionprocessed[1,"start"],end=regionprocessed[1,"end"])

namecoord<-paste("data/tempfilefst/",name,"_coordinates.txt",sep="")

write.table(coordinates,file=namecoord,quote=FALSE,row.names=FALSE,col.names = FALSE,sep="\t")


