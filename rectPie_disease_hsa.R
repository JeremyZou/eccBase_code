rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_disease/")
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggsunburst)
library(data.table)
library(ggnewscale)
library(gridExtra)
##################
rawF <- read.csv("../raw/hsa2.txt",sep = "\t",stringsAsFactors = F,header = T)
aa <- rawF[,c("Disease","Cell.Tissues","Identify_method")]
ad <- read.csv("dtcm",sep = "\t",header = F,stringsAsFactors = F)
df <- data.frame()
ti = 0
for (line in rownames(ad)){
  ti = ti+1
  print (ti)
  cc=paste0("row_",ti)
  print(cc)
  cc <- aa[grep(pattern = ad[line,]$V1, aa$Disease),]
  cc$Disease <- ad[line,]$V1
  cc <- cc[grep(pattern = ad[line,]$V2, cc$Cell.Tissues),]
  cc$Cell.Tissues <- ad[line,]$V2
  cc <- cc[which(cc$Identify_method == ad[line,]$V3),]
  cc$Identify_method <- ad[line,]$V3
  df <- rbind(df, cc)
}
dfd <- as.data.frame(table(df$Disease))
dfc <- as.data.frame(table(df$Cell.Tissues))
dfm <- as.data.frame(table(df$Identify_method))
df01 <- df
ti=0
for(kin in rownames(dfd)){
  ti=ti+1; print(ti)
  df01$Disease <- gsub(dfd$Var1[ti],dfd$Freq[ti],df01$Disease)
}
ti=0
for(kin in rownames(dfc)){
  ti=ti+1; print(ti)
  df01$Cell.Tissues <- gsub(dfc$Var1[ti],dfc$Freq[ti],df01$Cell.Tissues)
}
ti=0
for(kin in rownames(dfm)){
  ti=ti+1; print(ti)
  df01$Identify_method <- gsub(dfm$Var1[ti],dfm$Freq[ti],df01$Identify_method,
                               fixed = T)
}
df01$Identify_method <- gsub("74/ 3160","382",df01$Identify_method,fixed = T)
df_f<-cbind(df,df01)
colnames(df_f) <- c("disease","celltissue","imeth","couD","couC","couM")
df_f[grepl("Predict|Manually]",df_f$imeth),3] <- "AmpliconArchitect &| IGV"
write.table(df_f,"df_f",col.names = T,row.names = F,quote = F,sep = "\t")
df_f[which(as.numeric(df_f$couD) < 5000),c(1,2,3)] <- "other"
sum(as.numeric(unique(df_f[which(as.numeric(df_f$couD) < 5000),4])))

sum(as.numeric(unique(df_tissue[which(as.numeric(df_tissue$couC) > 5000),4])))
##################
#sum(as.numeric(unique(df_f[grepl(x = df_f[,1],pattern = "other"),][,4])))
#19872/750000
##################
Df_f <- as.data.table(df_f)
sb_0 <- Df_f[, .(xmin = 0.0, xmax = 1.2, ymin = 0, ymax = .N)]
sb_1 <- Df_f[, .N, by = c("disease")]
#排序方法一
objord <- c("other","Testicular germ cell tumor","Glioma","Stomach cancer",
            "Kidney cancer","Head and Neck squamous cell carcinoma",
            "Esophageal carcinoma","Liver cancer","Lung cancer","Breast cancer",
            "Colon adenocarcinoma","Cervical cancer","Lymphoma",
            "Prostate cancer","Ovarian cancer","Healthy")
sb_1$disease <- factor(sb_1$disease, levels = objord)
sb_1 <- sb_1[order(sb_1$disease),]
#排序方法二
#setkey(sb_1, N)
sb_1[, `:=`(xmin = 1.2, xmax = 2.5, ymin = cumsum(shift(N, fill = 0)), ymax = cumsum(N))]
sb_1$frac <- paste(round(sb_1$N/sum(sb_1$N),4)*100,"%")
#################
#start to plot
g_rec <-
  ggplot() +
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) +
  geom_rect(data = sb_0, fill = "white") +
  geom_rect(data = sb_1, mapping = aes(fill = disease), colour="white") +
  scale_fill_manual(breaks = objord,
                    values = rev(c(pal_d3(palette = "category20")(20)[1:16])),
                    labels = paste(c("Other","Testicular germ cell tumor","Glioma","Stomach cancer","Kidney cancer","Head and Neck squamous \n cell carcinoma","Esophageal carcinoma","Liver cancer","Lung cancer","Breast cancer","Colon adenocarcinoma","Cervical cancer","Lymphoma","Prostate cancer","Ovarian cancer","Healthy")," (",sb_1$frac,")")) +
#  geom_text_repel(data = sb_1, 
#                  aes(x=3.5,y=(ymax-ymin)/2,label=frac,hjust=0,
#                                   colour="black",family="serif"),size=4) +
  guides(fill = guide_legend(order = 1, nrow = 16, title = "Disease",reverse = T))
g_rec
g_sb_01 <-
  g_rec +
  coord_polar(theta = "y",direction = -1) +
  ggtitle(label = "hsa eccDNA distribution by disease") +
  theme_classic() +
  theme(legend.position = c(1.3,0.5),
        legend.text = element_text(family = "serif", size = 16),
        legend.title = element_text(family = "serif", size = 16),
        plot.margin = unit(c(0,10,0,0),"cm"),
        plot.title = element_text(family = "serif", size = 18,hjust = 0.5)) +
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank()) +
  annotate("text",label = "hsa (736775)\n---------------\nIn total", x=0, y=0.5,
           family="serif",size=6)
g_sb_01
#ggsave("rectPie_disease_hsa.jpg",dpi=800)
ggsave(g_sb_01,filename = "AI/rectPie_disease_hsa.pdf",units = "cm", width = 24, height = 15)











##################
##################
##################
#按disease数量取子集后，couD不需要重新计算，couC和couM需要
df_other <- df_f[which(as.numeric(df_f$couD)< 10000),]
df_oth2 <- df_other
liTc <- as.data.frame(table(df_other$celltissue))
liTm <- as.data.frame(table(df_other$imeth))
ti=0
for(kin in rownames(liTc)){
  ti=ti+1; print(ti)
  df_other$celltissue <- gsub(liTc$Var1[ti],liTc$Freq[ti],df_other$celltissue)
}
ti=0
for(kin in rownames(liTm)){
  ti=ti+1; print(ti)
  df_other$imeth <- gsub(liTm$Var1[ti],liTm$Freq[ti],df_other$imeth)
}
df_other$imeth <- gsub("74/ 1652","382",df_other$imeth,fixed = T)
df_bb<-cbind(df_oth2[,1:3],df_other[,c(4,2,3)])
df_bb$nameD <-"disss";df_bb$nameC <-"ccell";df_bb$nameM <-"metho"
colnames(df_bb) <- c("disease","celltissue","imeth","couD","couC",
                     "couM","nameD","nameC","nameM")
ti=0
for(line in rownames(df_bb)){
  ti=ti+1
  print(ti)
  df_bb$per_d[ti] <- round(as.numeric(df_bb$couD)[ti]/length(df_bb$couD),6)
  df_bb$per_c[ti] <- round(as.numeric(df_bb$couC)[ti]/length(df_bb$couD),6)
  df_bb$per_m[ti] <- round(as.numeric(df_bb$couM)[ti]/length(df_bb$couD),6)
}
df_bb<-df_bb[order(as.numeric(df_bb$couD),decreasing = T),]



