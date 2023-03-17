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
##################
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
  df01$Identify_method <- gsub(dfm$Var1[ti],dfm$Freq[ti],df01$Identify_method,fixed = T)
}
df01$Identify_method <- gsub("74/ 3160","382",df01$Identify_method,fixed = T)
df_f<-cbind(df,df01)
colnames(df_f) <- c("disease","celltissue","imeth","couD","couC","couM")
df_f[grepl("Predict|Manually]",df_f$imeth),3] <- "AmpliconArchitect &| IGV"
write.table(df_f,"df_f",col.names = T,row.names = F,quote = F,sep = "\t")
df_f[which(as.numeric(df_f$couD) < 10000),c(1,2,3)] <- "other"
#查看
length(table(df_f$disease))
table(df_f$imeth)
#作图
#Df_f <- as.data.table(df_f[grep("Ovarian", df_f$disease),])
Df_f <- as.data.table(df_f)
sb_0 <- Df_f[, .(xmin = 0.0, xmax = 0.8, ymin = 0, ymax = .N)]   # Initial Data Set: just count the records
sb_1 <- Df_f[, .N, by = c("disease")]                            # Counts by category 1
sb_2 <- Df_f[, .N, by = c("disease", "celltissue")]              # Counts by category 1 and 2
sb_3 <- Df_f[, .N, by = c("disease", "celltissue", "imeth")]     # Counts by category 1, 2, and 3
write.table(sb_3,"sb_3",sep = "\t",quote = F,row.names = F,col.names = F)





# insure that the data sets are ordered by category 1, 2, then 3
setkey(sb_1, disease)
setkey(sb_2, disease, celltissue)
setkey(sb_3, disease, celltissue, imeth)

sb_1[, `:=`(xmin = 0.8, xmax = 2.5, ymin = cumsum(shift(N, fill = 0)), ymax = cumsum(N))]
sb_2[, `:=`(xmin = 2.5, xmax = 4.2, ymin = cumsum(shift(N, fill = 0)), ymax = cumsum(N))]
sb_3[, `:=`(xmin = 4.2, xmax = 4.6, ymin = cumsum(shift(N, fill = 0)), ymax = cumsum(N))]
g_rec <-
  ggplot() +
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) +
  geom_rect(data = sb_0, fill = "white") +
  
  # add the layers for the cut
  geom_rect(data = sb_1, mapping = aes(fill = disease), colour="white") +
  #scale_fill_brewer(type = "qual", palette = 1) +
  #scale_fill_d3("category20")+
  #scale_fill_igv()+
  scale_fill_simpsons() +
  
  # allow for a new fill scale, add the layer for the color by cut
  new_scale("fill") +
  geom_rect(data = sb_2, colour ="white", mapping = aes(fill = celltissue)) +
  # scale_fill_brewer(type = "qual", palette = 2) +
  scale_fill_d3("category20") +
  # allow for a new fill scale, add the layer for the clarity by color by cut
  new_scale("fill") +
  geom_rect(data = sb_3, mapping = aes(fill = imeth),colour="white") +
  #scale_fill_brewer(palette = 1) +
  #scale_fill_d3("category20") +
  #scale_fill_simpsons() +
  #scale_fill_manual(values=c("mediumpurple","indianred","gray60","grey","yellow")) +
  
  theme(legend.position = "bottom") +
  guides(fill         = guide_legend(order = 5, nrow = 5),
         fill_new     = guide_legend(order = 4, nrow = 5),
         fill_new_new = guide_legend(order = 3, nrow = 5))
#guide_legend()
#g_rec
# build the sunburst plot
g_sb <-
  g_rec +
  coord_polar(theta = "y",direction = -1) +
  theme_classic() +
  theme(legend.position = "bottom",
        axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank())
g_sb
ggsave("test.jpg",dpi=800)
View(sb_3)











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
##################
##################
##################




