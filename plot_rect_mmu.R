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
rawF <- read.csv("./dtcm_mmu",sep = "\t",stringsAsFactors = F,header = F)
#rawF <- rawF[,c(4,1,2,3)]
colnames(rawF) <- c("disease","celltissue","imeth","NN")
Df_f <- as.data.table(rawF)
sum(Df_f$NN)

#作图
#Df_f <- as.data.table(df_f[grep("Ovarian", df_f$disease),])

sb_0 <- Df_f[, .(xmin = 0.0, xmax = 0.8, ymin = 0, ymax = .N)]
sb_1 <- data.frame(disease="Healthy",NN=sum(Df_f$NN)); sb_1 <- as.data.table(sb_1)
sb_2 <- Df_f[, c("disease", "celltissue","NN")]
sb_3 <- Df_f




# insure that the data sets are ordered by category 1, 2, then 3
setkey(sb_1, disease)
setkey(sb_2, disease, celltissue)
setkey(sb_3, disease, celltissue, imeth)

sb_1[, `:=`(xmin = 0.8, xmax = 2.5, ymin = cumsum(shift(NN, fill = 0)), ymax = cumsum(NN))]
sb_2[, `:=`(xmin = 2.5, xmax = 4.2, ymin = cumsum(shift(NN, fill = 0)), ymax = cumsum(NN))]
sb_3[, `:=`(xmin = 4.2, xmax = 4.6, ymin = cumsum(shift(NN, fill = 0)), ymax = cumsum(NN))]
g_rec <-
  ggplot() +
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) +
  geom_rect(data = sb_0, fill = "white") +
  
  # add the layers for the cut
  geom_rect(data = sb_1, mapping = aes(fill = disease)) +
  #scale_fill_brewer(type = "qual", palette = 1) +
  #scale_fill_d3("category20")+
  #scale_fill_igv()+
  scale_fill_manual(values = "red") +
  
  # allow for a new fill scale, add the layer for the color by cut
  new_scale("fill") +
  geom_rect(data = sb_2, colour ="white", mapping = aes(fill = celltissue)) +
  # scale_fill_brewer(type = "qual", palette = 2) +
  #scale_fill_d3("category20") +
  scale_fill_igv() +
  # allow for a new fill scale, add the layer for the clarity by color by cut
  new_scale("fill") +
  geom_rect(data = sb_3, mapping = aes(fill = imeth)) +
  #scale_fill_brewer(palette = 1) +
  #scale_fill_d3("category20") +
  #scale_fill_simpsons() +
  #scale_fill_manual(values=c("mediumpurple","indianred","gray60","grey","yellow")) +
  
  theme(legend.position = "bottom") +
  guides(fill         = guide_legend(order = 5, nrow = 5),
         fill_new     = guide_legend(order = 4, nrow = 5),
         fill_new_new = guide_legend(order = 3, nrow = 5))

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




