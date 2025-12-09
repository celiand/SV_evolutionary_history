#### To use this script, first save PCA data from the summarize_run lostruct scripts. Then eventually associate individual with their population
#### The input data is then: Row name = individual ID | PC1 | PC2 | population | window (window being the window number from the lostruct package that cuts the chromosomes into several windows)

data_pca<-read.table(file="yourdata")

#### this function try to make the best cluster for different number of k
cluster_windows<-function(window){
  subdata<-data_pca[data_pca$window==window,c(1,2)]
  values<-c()
  for(i in 1:10){
    newk<-c()
    for(y in 1:100){
      test<-kmeans(subdata,centers=i)
      newk<-c(newk,test$tot.withinss)
    }
    values<-c(values,mean(newk))
  }
  datakcluster<-data.frame(x=1:10,y=values)
  return(datakcluster)
}

### list all unique windows in a specific chromosome for example
allwind<-unique(data_pca$window)

### choose a window to investigate
window<-allwind[10]

### make cluster file and have a look at the kmeans
datakcluster<-cluster_windows(window = window)

ggplot(datakcluster,aes(x=x,y=y))+geom_line()+scale_x_continuous(breaks=seq(1,10,1))

### select a number of cluster based on the kmeans
nbclusters<-6



### plot the groups on the PCA
datatested<-data_pca[data_pca$window==window,c(1,2)]
test<-kmeans(datatested,centers=nbclusters)

datatested$cluster<-test$cluster

ggplot(datatested, aes(x = PC1, y = PC2,color=as.factor(cluster))) +
  geom_point() +theme_bw(20)


### do several tries. Sometimes the kmeans have trouble to group individuals, especially when they are very grouped on 1 axis but very spread on the other

### when sometimes seems correct, you can save the data (which individual belongs to which cluster)
for( i in 1:length(datatested$cluster)){
  name<-unlist(strsplit(rownamescorrect[i],split="_"))
  startname<-name[3]
  for( y in 4:length(name)){
    if(name[y]==startname){
      samplename<-paste(name[3:(y-1)],collapse ="_")
      datatested[i,"ind"]<-samplename
      break
    }
  }
}


datatested$start<-regionspca[window,"start"]-1000000
datatested$end<-regionspca[window,"end"]+1000000
datatested$chrom<-regionspca[window,"chrom"]
datatested$window<-window
datatested$code<-paste(unique(datatested$chrom),unique(datatested$start),unique(datatested$end),nbclusters,sep="_")

### import previous data saved and add the new window
previousdata<-read.table(file="full_info_windows_PCA.txt",header=TRUE)


newdatatested<-rbind.data.frame(previousdata,datatested,stringsAsFactors = FALSE)

write.table(newdatatested,file="full_info_windows_PCA.txt",row.names = FALSE)


#### at the end, the table generated contain region coordinate individual names and their group in each regions of interest
#### it's then possible to make a FST between each group for each region