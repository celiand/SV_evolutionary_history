allfile<-list.files("data/fst_files_cluster/")

library(googledrive)
library(ggplot2)
##make a loop to plot all fst files and register them into google drive "group_fst" (for instance)

for( files in allfile){
  filename<-paste("data/fst_files_cluster/",files,sep="")
  data<-read.table(file=filename,header=TRUE)
  plotwind<-ggplot(data,aes(x=POS,y=WEIR_AND_COCKERHAM_FST))+geom_point(alpha=0.6,size=3)+theme_bw(18)
  fname_export<-paste(files,".png",sep="_")
  ggsave(filename=fname_export,plot=plotwind,device="png",path="data/analysis/fst_pca",width=20,height=16,units="in")
  pathlocal<-paste("data/analysis/fst_pca",fname_export,sep="/")
  drive_upload(media=pathlocal,name=fname_export,path="group_fst",overwrite=TRUE)
}