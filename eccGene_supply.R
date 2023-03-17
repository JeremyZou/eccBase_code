#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_gene/")
#install.packages("UpSetR")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","ggupset","tidyverse",
         "UpSetR")
for (i in pks){library(i, character.only = T)}
rm(i,pks)

#开始整理数据：智人的全部基因类型count字段的upsetR数据
ha <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/hsa2.txt", sep="\t", stringsAsFactors = F, header = T)
hb <- cbind(ha[, grepl("count",colnames(ha))], ha[, c(1, 5, 6, 48)])
hG <- hb[,1:5]
hG$nnam <- as.numeric(rownames(hG))
ti = 0
for (i in seq(length(colnames(hG))-1)) {
  ti = ti + 1; print(ti)
  hG[,i] <- apply(hG, 1, function(x){
  ifelse(x[i]==0, x[i] <- 0, x[i] <- x[6])
})
}
h_listInput <- list(Protein_coding    = hG[which(hG[,2] >0),2], 
                    ncRNA      = hG[which(hG[,3] >0),3],  
                    lncRNA     = hG[which(hG[,4] >0),4], 
                    pseudogene = hG[which(hG[,5] >0),5])

#开始整理数据：小鼠的全部基因类型count字段的upsetR数据
ma <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/mmu.txt", sep="\t", stringsAsFactors = F, header = T)
mb <- cbind(ma[, grepl("count", colnames(ma))], ma[,c(1,5,6,48)])
mG <- mb[,1:5]
mG$nnam <- as.numeric(rownames(mG))
ti = 0
for (i in seq(length(colnames(mG))-1)) {
  ti = ti + 1; print(ti)
  mG[,i] <- apply(mG, 1, function(x){
    ifelse(x[i]==0, x[i] <- 0, x[i] <- x[6])
  })
}
m_listInput <- list(Protein_coding    = mG[which(mG[,2] >0),2], 
                    ncRNA      = mG[which(mG[,3] >0),3],  
                    lncRNA     = mG[which(mG[,4] >0),4], 
                    pseudogene = mG[which(mG[,5] >0),5])


#开始绘图：智人的全部基因类型count字段的upsetR plot
jpeg(filename = "./supply/geneSupply_hsa_blueBlue.jpg", res = 800,units = "cm", width = 28, height = 12)
pdf(file = "./supply/geneSupply_hsa_blueBlue.pdf", width = 14, height = 5)
upset(fromList(h_listInput), mb.ratio = c(0.75, 0.25), order.by = "freq",
      main.bar.color = "royalblue4", matrix.color = "royalblue4",text.scale = 1.2,
      mainbar.y.label = "Number of eccDNAs",sets.x.label = "Number of eccDNAs",
      point.size = 2.5, line.size = 1, sets.bar.color = "royalblue4",
      set_size.show = T,set_size.scale_max = 736775)

dev.off()

#开始绘图：小鼠的全部基因类型count字段的upsetR plot
jpeg(filename = "./supply/geneSupply_mmu_blueBlue.jpg", res = 800,units = "cm", width = 28, height = 12)
pdf(file = "./supply/geneSupply_mmu_blueBlue.pdf", width = 14, height = 5)
upset(fromList(m_listInput), mb.ratio = c(0.75, 0.25), order.by = "freq",
      main.bar.color = "royalblue4", matrix.color = "royalblue4",text.scale = 1.2,
      mainbar.y.label = "Number of eccDNAs",sets.x.label = "Number of eccDNAs",
      point.size = 2.5, line.size = 1, sets.bar.color = "royalblue4",
      set_size.show = T,set_size.scale_max = 400000)

dev.off()











