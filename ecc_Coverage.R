###测试测试测试
rm(list=ls())
setwd("/disk5/shy/5eccdna/ecc_density/")
library(karyoploteR)
###载入数据
aa <- read.csv("/disk5/shy/5eccdna/raw/hsa.txt",sep = "\t",
               stringsAsFactors = F, header = TRUE)
aa_d <- aa[grep(pattern = "Healthy",aa[,"Disease"]),c("Coordinate.hg19.","eccDNA_ID","Disease")]
aa_d <- aa_d[which(!(is.na(aa_d$Coordinate.hg19.))),]
aa_d[,4:6] <- as.data.frame(str_split_fixed(aa_d$Coordinate.hg19., ":|-", 3))
aa_d <- aa_d[,c(4,5,6,2)]
#aa_d <- read.table(file = "/disk5/shy/1shyJOB/other/healthy.bed",sep = "\t",stringsAsFactors = F)
colnames(aa_d) <- c("chr","start","end","info")
aa_G <- makeGRangesFromDataFrame(aa_d)
###开始作图
#jpeg("kpPlotCoverage_Healthy.jpg",res=400, width = 2000, height = 3000)
jpeg("kpPlotCoverage_Cancer.jpg",res=400, width = 2000, height = 3000)
kp <- plotKaryotype(genome="hg19", plot.type=1, cex=1)
kpAxis(kp, cex=0.01,numticks=2, labels=NA)
kp <- kpPlotCoverage(kp, data=aa_G,show.0.cov = F)
#kp <- kpPlotDensity(kp, data=aa_G)
#bw_f <- "minsorthg38.bw"
#bw_t <- "h3.bw"
#kp <- kpPlotBigWig(kp, data=bw_t, show.0.cov=FALSE, col="hotpink3")
#kp
dev.off()







