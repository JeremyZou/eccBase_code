#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_predict/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling")
for (i in pks){library(i, character.only = T)}
rm(i,pks)

#开始整理数据：导入数据、拆分组别、以及标准化
ha   <- read.csv(file = "../raw/hsa2.txt", sep = "\t", stringsAsFactors = F, header = T)
ha_h <- ha[grep(pattern = "Healthy", ha[,"Disease"]), c("eccDNA_ID", "Cell.Tissues", "Length", colnames(ha)[grepl("frac",colnames(ha))])]
ha_h <- na.omit(ha_h)
for (i in seq(3,length(colnames(ha_h)))){
  print(i)
  mm=mean(ha_h[,i]);mx=max(ha_h[,i]);mi=min(ha_h[,i]);sd=sd(ha_h[,i])
  ha_h[,i] <- sapply(ha_h[,i],function(x){(x-mm)/sd})
}
write.table(ha_h,file = "ha_h.csv", quote = F, sep = ",",col.names = T,row.names = F)

ha_c <- ha[grep(pattern = "[^Healthy]", ha[,"Disease"]), c("eccDNA_ID", "Disease", "Length", colnames(ha)[grepl("frac",colnames(ha))])]
ha_c <- na.omit(ha_c)
for (i in seq(3,length(colnames(ha_c)))){
  print(i)
  mm=mean(ha_c[,i]);mx=max(ha_c[,i]);mi=min(ha_c[,i]);sd=sd(ha_h[,i])
  ha_c[,i] <- sapply(ha_c[,i],function(x){(x-mm)/sd})
}
write.table(ha_c,file = "ha_c.csv", quote = F, sep = ",",col.names = T,row.names = F)

#开始整理数据：联合Cancer、Healthy
ha_h2 <- ha_h; ha_h2$eccDNA_ID <- 0
ha_c2 <- ha_c; ha_c2$eccDNA_ID <- 1
ha_u <- rbind(ha_c2,ha_h2)
write.table(ha_u,file = "ha_uALL.csv",sep = ",",quote = F,col.names = T,row.names = F)

#开始整理数据：随机抽取Cancer、Healthy各自10万数据
set.seed(111)
ha_h3 <- ha_h2[base::sample(rownames(ha_h2),100000),]
ha_c3 <- ha_c2[base::sample(rownames(ha_c2),100000),]
ha_u2 <- rbind(ha_c3,ha_h3)
write.table(ha_u2,file = "ha_u100000.csv",sep = ",",quote = F,col.names = T,row.names = F)

#开始整理数据：合并Disease的Glioma、以及随机抽取Healthy的Blood等量数据
set.seed(111)
ha_Glioma <- ha_c2[grep(pattern = "Glioma", ha_c2[,"Disease"]), c(1,seq(3,length(colnames(ha_c2))))]
ha_hG     <- ha_h2[base::sample(rownames(ha_h2),length(ha_Glioma[,1])),]
ha_u3     <- rbind(ha_Glioma,ha_hG[,-c(2)])
write.table(ha_u3,file = "ha_Glioma.csv",sep = ",",quote = F,row.names = F,col.names = T)


#开始整理数据：合并Disease的Glioma、以及Healthy的Blood
ha_Glioma <- ha_c2[grep(pattern = "Glioma", ha_c2[,"Disease"]), c(1,seq(3,length(colnames(ha_c2))))]
ha_Blood  <- ha_h[grep(pattern = "Blood", x = ha_h$Cell.Tissues),c(1,seq(3,18))]
ha_Blood$eccDNA_ID <- 0
ha_u4     <- rbind(ha_Glioma,ha_Blood)
write.table(ha_u4,file = "ha_GliomaBlood.csv",sep = ",",quote = F,row.names = F,col.names = T)


#开始整理数据：合并Disease的Stomach cancer、以及Healthy的Blood
ha_Stomach <- ha_c2[grep(pattern = "Stomach", ha_c2[,"Disease"]), c(1,seq(3,length(colnames(ha_c2))))]
ha_Blood   <- ha_h[grep(pattern = "Blood", x = ha_h$Cell.Tissues),c(1,seq(3,18))]
ha_Blood$eccDNA_ID <- 0
ha_u5      <- rbind(ha_Stomach,ha_Blood)
write.table(ha_u5,file = "ha_StomachBlood.csv",sep = ",",quote = F,row.names = F,col.names = T)
set.seed(111)
ha_u5Test  <- rbind(ha_Stomach[base::sample(rownames(ha_Stomach), replace = F, size = 1022),],
                    ha_Blood[base::sample(rownames(ha_Blood), replace = F, size = 897),])
ha_u5Train <- ha_u5[which(!(rownames(ha_u5) %in% rownames(ha_u5Test))),]
write.table(ha_u5Test, file = "ha_StomachBlood_test.csv", sep = ",",quote = F,row.names = F,col.names = T)
write.table(ha_u5Train,file = "ha_StomachBlood_train.csv",sep = ",",quote = F,row.names = F,col.names = T)


#开始整理数据：预先标准化
ha_t <- ha[, c("eccDNA_ID", "Disease","Cell.Tissues", "Length", colnames(ha)[grepl("frac",colnames(ha))])]
for (i in seq(4,length(colnames(ha_t)))){
  print(i)
  mm=mean(ha_t[,i]);mx=max(ha_t[,i]);mi=min(ha_t[,i]);sd=sd(ha_t[,i])
  ha_t[,i] <- sapply(ha_t[,i],function(x){(x-mm)/sd})
}

#开始整理数据：预先标准化 + 合并Disease的Glioma、以及Healthy的Blood
ha_Glioma <- ha_t[grep(pattern = "Glioma", ha_t[,"Disease"]), c(1,seq(4,length(colnames(ha_t))))]
ha_Blood  <- ha_t[grep(pattern = "Blood",  ha_t$Cell.Tissues),c(1,seq(4,length(colnames(ha_t))))]
ha_Blood$eccDNA_ID <- 0; ha_Glioma$eccDNA_ID <- 1
ha_u4     <- rbind(ha_Glioma, ha_Blood)
write.table(ha_u4,file = "ha_GliomaBlood.csv",sep = ",",quote = F,row.names = F,col.names = T)

#开始整理数据：预先标准化 + 合并Disease的Stomach cancer、以及Healthy的Blood
ha_Stomach <- ha_t[grep(pattern = "Stomach", ha_t[,"Disease"]), c(1,seq(4,length(colnames(ha_t))))]
ha_Blood   <- ha_t[grep(pattern = "Blood",   ha_t$Cell.Tissues),c(1,seq(4,length(colnames(ha_t))))]
ha_Blood$eccDNA_ID <- 0; ha_Stomach$eccDNA_ID <- 1
ha_u4     <- rbind(ha_Stomach, ha_Blood)
write.table(ha_u4,file = "ha_StomachBlood.csv", sep = ",", quote = F, row.names = F, col.names = T)

#开始整理数据：预先标准化 + 合并Cell.Tissues的PC-3、以及Healthy的Blood
ha_PC    <- ha_t[grep(pattern = "PC-3", ha_t[,"Cell.Tissues"]), c(1,seq(4,length(colnames(ha_t))))]
ha_Blood <- ha_t[grep(pattern = "Blood",   ha_t$Cell.Tissues),c(1,seq(4,length(colnames(ha_t))))]
ha_Blood$eccDNA_ID <- 0; ha_PC$eccDNA_ID <- 1
ha_u4    <- rbind(ha_PC, ha_Blood)
write.table(ha_u4,file = "ha_PC3Blood.csv", sep = ",", quote = F, row.names = F, col.names = T)

write.table(ha_u4,file = "ha_StomachBlood.csv", sep = ",", quote = F, row.names = F, col.names = T)


