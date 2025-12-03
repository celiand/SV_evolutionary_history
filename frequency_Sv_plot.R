## ----  Start to plot Lostruct_Sv frequency ---- 

regionspca<-read.table(file="Data/interesting_region_bis.txt",header=TRUE)



all_data<-c()
for (i in 1:nrow(regionspca)){
  
  filename_base<-paste(regionspca[i,"CHROM"],regionspca[i,"Manual_start"],regionspca[i,"Manual_end"],sep="_")
  
  filename<-paste("Data/sprgi_gt_pca_files/",filename_base,".txt",sep="")  
  
  filetest<-read.table(file=filename,header=TRUE)
  filetest$genonb<-filetest$GT
  filetest[filetest$genonb=="NI/NI","genonb"]<-0
  filetest[filetest$genonb=="NI/I","genonb"]<-1
  filetest[filetest$genonb=="I/I","genonb"]<-2
  
  filetest$genonb<-as.numeric(as.character(filetest$genonb))
  
  
  freq_general<-sum(filetest$genonb)/740
  freq_FA<-sum(filetest[filetest$pop=="Farmed_American","genonb"])/160
  freq_FE<-sum(filetest[filetest$pop=="Farmed_European","genonb"])/224
  freq_WA<-sum(filetest[filetest$pop=="Wild_American","genonb"])/160
  freq_WE<-sum(filetest[filetest$pop=="Wild_European","genonb"])/196
  
  
  
  
  data<-data.frame(Pop=c("All","Farmed_American","Farmed_European","Wild_American","Wild_European"),freq=c(freq_general,freq_FA,freq_FE,freq_WA,freq_WE),continent=c("All_data","North_America","Europe","North_America","Europe"),type=c("All_data","Farmed","Farmed","Wild","Wild"))
  data$region<-regionspca[i,"newcode"]
  
  all_data<-rbind.data.frame(all_data,data,stringsAsFactors = FALSE)
  
}


### Plot frequency within continent
library(ggplot2)
freqallAM<-all_data[all_data$continent=="North_America",]
freqallEU<-all_data[all_data$continent=="Europe",]
freqall<-data.frame(regions=freqallAM$region,freqAM=freqallAM$freq,freqEU=freqallEU$freq,type=freqallAM$type)
freqall$color<-"black"
freqall[5:6,"color"]<-"firebrick1" ###color the outlier


ggplot(freqall,aes(x=freqAM,y=freqEU,color=color))+theme_bw(20)+geom_point(size=3)+geom_smooth(method = "lm", se = TRUE)+facet_wrap(~type)+scale_color_identity()

### Plot frequency for each domestication status
freqallWild<-all_data[all_data$type=="Wild",]
freqallFarm<-all_data[all_data$type=="Farmed",]
freqall_bis<-data.frame(regions=freqallWild$region,freqWild=freqallWild$freq,freqFarmed=freqallFarm$freq,continent=freqallWild$continent)
freqall_bis$color<-"black"
freqall_bis[5:6,"color"]<-"firebrick1" ###color the outlier

ggplot(freqall_bis,aes(x=freqWild,y=freqFarmed,color=color))+geom_point(size=3)+theme_bw(20)+facet_wrap(~continent)+geom_smooth(method = "lm", se = TRUE)+scale_color_identity()

## make statistical test for the correlation with PCA SV only

## adding the type or continent as a fixed effect
#model <- lm(freqAM ~ freqEU + type, data = freqall)
#summary(model)

#model2 <- lm(freqWild ~ freqFarmed + continent, data = freqall_bis)
#summary(model2)


### within each group test
table_wild<-freqall[freqall$type=="Wild",]
model <- lm(freqAM ~ freqEU, data = table_wild)
summary(model)

table_farmed<-freqall[freqall$type=="Farmed",]
model2 <- lm(freqAM ~ freqEU, data = table_farmed)
summary(model2)



table_europe<-freqall_bis[freqall_bis$continent=="Europe",]
model3 <- lm(freqWild ~ freqFarmed, data = table_europe)
summary(model3)

table_america<-freqall_bis[freqall_bis$continent=="North_America",]
model4 <- lm(freqWild ~ freqFarmed, data = table_america)
summary(model4)


## ---- add bayestyper_SV as background data ---- 


FE_SV<-read.table(file="Data/bayestyper_SV_info_table_Farmed_European.txt",header=TRUE)
FE_SV$continent<-"Europe"
FE_SV$type<-"Farmed"

WE_SV<-read.table(file="Data/bayestyper_SV_info_table_Wild_European_nooutlier.txt",header=TRUE)
WE_SV$continent<-"Europe"
WE_SV$type<-"Wild"

WA_SV<-read.table(file="Data/bayestyper_SV_info_table_Wild_American.txt",header=TRUE)
WA_SV$continent<-"American"
WA_SV$type<-"Wild"

FA_SV<-read.table(file="Data/bayestyper_SV_info_table_Farmed_American.txt",header=TRUE)
FA_SV$continent<-"American"
FA_SV$type<-"Farmed"

all_SV_BT<-rbind.data.frame(FE_SV,WE_SV,WA_SV,FA_SV,stringsAsFactors = FALSE)


library(ggplot2)
freqall_BTAM<-all_SV_BT[all_SV_BT$continent=="American",]
freqall_BTEU<-all_SV_BT[all_SV_BT$continent=="Europe",]
freqall_BT<-data.frame(freqAM=freqall_BTAM$FREQ,freqEU=freqall_BTEU$FREQ,type=freqall_BTAM$type)


model_BT <- lm(freqAM ~ freqEU + type, data = freqall_BT)
summary(model_BT)

freqall_BTWild<-all_SV_BT[all_SV_BT$type=="Wild",]
freqall_BTFarm<-all_SV_BT[all_SV_BT$type=="Farmed",]
freqall_BT_bis<-data.frame(freqWild=freqall_BTWild$FREQ,freqFarmed=freqall_BTFarm$FREQ,continent=freqall_BTWild$continent)





#### combined plots

ggplot(freqall_bis[freqall_bis$continent=="Europe",],aes(x=freqWild,y=freqFarmed))+theme_bw(20)+
  stat_density2d(
    aes(x = freqWild, y = freqFarmed, fill = ..density..^0.25), 
    data = freqall_BT_bis[freqall_BT_bis$continent == "Europe", ], 
    geom = "tile", 
    contour = FALSE, 
    n = 150,
    alpha = 0.5  
  ) +
  scale_fill_viridis_c(option = "G", begin = 0.1, end = 0.9, direction=-1) +geom_point(size=2,aes(color = color))+ geom_smooth(method = "lm", se = FALSE, color = "dodgerblue4")+
  scale_color_identity()+ggtitle("Freq Wild vs Domesticated - Europe")


ggplot(freqall_bis[freqall_bis$continent=="North_America",],aes(x=freqWild,y=freqFarmed))+theme_bw(20)+
  stat_density2d(
    aes(x = freqWild, y = freqFarmed, fill = ..density..^0.25), 
    data = freqall_BT_bis[freqall_BT_bis$continent == "American", ], 
    geom = "tile", 
    contour = FALSE, 
    n = 150,
    alpha = 0.5  
  ) +
  scale_fill_viridis_c(option = "G", begin = 0.1, end = 0.9, direction=-1) +geom_point(size=2,aes(color = color))+ geom_smooth(method = "lm", se = FALSE, color = "dodgerblue4")+
  scale_color_identity()+ggtitle("Freq Wild vs Domesticated - America")



ggplot(freqall[freqall$type=="Wild",],aes(x=freqEU,y=freqAM))+theme_bw(20)+
  stat_density2d(
    aes(x = freqEU, y = freqAM, fill = ..density..^0.25), 
    data = freqall_BT[freqall_BT$type=="Wild", ], 
    geom = "tile", 
    contour = FALSE, 
    n = 150,
    alpha = 0.5  
  ) +
  scale_fill_viridis_c(option = "G", begin = 0.1, end = 0.9, direction=-1) +geom_point(size=2,aes(color = color))+ geom_smooth(method = "lm", se = FALSE, color = "dodgerblue4")+
  scale_color_identity()+ggtitle("Freq Europe vs America - Wild")


ggplot(freqall[freqall$type=="Farmed",],aes(x=freqEU,y=freqAM))+theme_bw(20)+
  stat_density2d(
    aes(x = freqEU, y = freqAM, fill = ..density..^0.25), 
    data = freqall_BT[freqall_BT$type=="Farmed", ], 
    geom = "tile", 
    contour = FALSE, 
    n = 150,
    alpha = 0.5  
  ) +
  scale_fill_viridis_c(option = "G", begin = 0.1, end = 0.9, direction=-1) +geom_point(size=2,aes(color = color))+ geom_smooth(method = "lm", se = FALSE, color = "dodgerblue4")+
  scale_color_identity()+ggtitle("Freq Europe vs America - Domesticated")








