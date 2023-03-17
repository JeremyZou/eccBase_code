rm(list=ls())
aa <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/hsa.txt",
                 sep="\t",stringsAsFactors = F,header = T)
peakHealthy <- aa[grep(pattern = "Healthy",aa$Disease),c("Coordinate.hg19.","eccDNA_ID","Disease")]
#peakHealthy <- aa[grep(pattern = "Predict",aa$Identify_method),c("Coordinate.hg19.","eccDNA_ID","Identify_method")]
peak_aa <- na.omit(peakHealthy)
peak_aa[,4:6] <- as.data.frame(str_split_fixed(peak_aa$Coordinate.hg19.,":|-",3))
peak_aa <- peak_aa[,c(4:6,2)]
colnames(peak_aa) <- c("chrom","start","end","eccDNA_ID")
peak_aa <- makeGRangesFromDataFrame(peak_aa)
peakCancer <- aa[grep(pattern = "[^Healthy]",aa$Disease),c("Coordinate.hg19.","eccDNA_ID","Disease")]
peak_bb <- na.omit(peakCancer)
peak_bb[,4:6] <- as.data.frame(str_split_fixed(peak_bb$Coordinate.hg19.,":|-",3))
peak_bb <- peak_bb[,c(4:6,2)]
colnames(peak_bb) <- c("chrom","start","end","eccDNA_ID")
peak_bb <- makeGRangesFromDataFrame(peak_bb)
peak01 <- peak_aa; peak02 <- peak_bb
peaks <- list(peak1=peak01,peak2=peak02)
#peak_t <- as.data.frame(as.numeric(as.character(peak_aa[,2])))
#peak01 <- read.csv("/disk5/shy/1shyJOB/7eccdna/ecc_density/fileBigWig/healthy.bed",sep="\t",header = F,stringsAsFactors = F)
#peak02 <- read.csv("/disk5/shy/1shyJOB/7eccdna/ecc_density/fileBigWig/healthy.bed",sep="\t",header = F,stringsAsFactors = F)
#作图
require(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#peak01 <- GRanges(seqnames=Rle(peak_aa[,1]),ranges=IRanges(peak_aa[,2], peak_aa[,3]), strand=rep(c("*"), nrow(peak_aa)))
#peak02 <- GRanges(seqnames=Rle(peak_bb[,1]),ranges=IRanges(peak_bb[,2], peak_bb[,3]), strand=rep(c("*"), nrow(peak_bb)))
txdb=TxDb.Hsapiens.UCSC.hg19.knownGene
peakAnnoList <- lapply(peaks, annotatePeak, TxDb = txdb, tssRegion = c(-1000, 1000))
write.csv(as.data.frame(as.GRanges(peakAnnoList[[1]])), 
          "/disk5/shy/1shyJOB/7eccdna/ecc_region/hsa_Healthy_chipseeker_anno.csv")
#pdf("/disk5/shy/1shyJOB/7eccdna/ecc_region/plotAnnoBar.pdf",width = 12,height = 5)
jpeg("/disk5/shy/1shyJOB/7eccdna/ecc_region/plotAnnoBar.jpg",width = 24,height = 8,res = 600,units = "cm")
#peakAnno <- annotatePeak(peak01, tssRegion=c(-1000, 1000),TxDb=txdb, annoDb="org.Hs.eg.db")
#pos_anno=as.data.frame(peakAnno)
#plotAnnoPie(peakAnno,legend.position = "rightside")
#plotAnnoPie(peakAnnoList[[1]])
#vennpie(peakAnno,r=0.25)
#vennpie(peakAnnoList[[1]])
#vennpie(peakAnnoList[[1]])
#plotDistToTSS(peakAnno)
plotAnnoBar(peakAnnoList)
#plotAnnoBar(peakAnno)
#plotDistToTSS(peakAnnoList)
dev.off()


