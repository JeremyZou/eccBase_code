#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling","IRanges",
         "GenomicRanges","rtracklayer","data.table","patchwork")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入TE的数据
dh <- c("LTR","LINE","SINE")
repeat_infor = "/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/mm10_rmsk_TE.gtf"
rpeat <- read.delim(repeat_infor,header=F,stringsAsFactors = F)
colnames(rpeat) <- c("chromosome","a","b","start","end","vaue","strand","oth","feature")
rpeat$len <- rpeat$end - rpeat$start
rpeatx <- rpeat[which(grepl(dh[3], rpeat$feature)), c(1, 4, 5, 9, 10)]
rpeatx <- makeGRangesFromDataFrame(rpeatx, keep.extra.columns = T)
#开始整理数据：与eccDNA进行findOverlap
f_all <- list.files(path = "./tempTissues_mmu", pattern = ".txt", full.names = T)
aa = ".*"
#aa=c("ERV1","ERVK","ERVL","ERVL-MaLR","Gypsy","LTR")
repFam_frac <- data.frame()
for (i in seq(length(list.files(path = "./tempTissues_mmu", pattern = ".txt", full.names = T)))) {
  print(i)
  temp  <- read.csv(f_all[i], header = F, sep = "\t", stringsAsFactors = F)
  temp  <- temp[,c(1,2,4)]
  temp  <- na.omit(temp)
  temp  <- separate(data = temp, col = V2, into = c("chromosome", "start","end"), sep = "[-|:]")
  temp[,c(3)] <- as.numeric(temp[,c(3)])
  temp[,c(4)] <- as.numeric(temp[,c(4)])
  temp  <- makeGRangesFromDataFrame(temp, keep.extra.columns = T)
  for (k in aa) {
    print(k)
    inter <- GenomicRanges::intersect(temp, 
                                      rpeatx[which(grepl(pattern = paste0("family_id ",k,";"), x = rpeatx$feature))], 
                                      ignore.strand = T)
    inter <- as.data.frame(inter)
    inter$feature <- inter$width - 1
    inter <- makeGRangesFromDataFrame(inter, keep.extra.columns = T)
    EDG   <- GenomicRanges::findOverlaps(temp, inter, ignore.strand = T)
    qwe   <- data.frame(eccN <- temp[queryHits(EDG)]$V1,
                        eccL <- temp[queryHits(EDG)]$V4,
                        ltrL <- inter[subjectHits(EDG)]$feature)
    colnames(qwe) <- c("eccN","eccL","ltrL")
    qwe$frac <- round(qwe$ltrL/qwe$eccL, 3)
    qwe <- as.data.table(qwe)
    qwe <- qwe[,.(frac.sum = sum(frac)),by=eccN]
    data <- data.frame(mark <- rep(str_split(basename(f_all[i]), pattern = "\\.")[[1]][1],
                                   length(temp$V4)),
                       frac <- c(qwe$frac.sum,rep(0,length(temp$V4)-length(qwe$frac.sum))))
    colnames(data) <- c("mark","frac")
    repFam_frac <- rbind(repFam_frac,data)
  }
}

kruskal.test(data=repFam_frac, frac~mark)
ccb <- with(data=repFam_frac, pairwise.wilcox.test(x=frac,g=mark,p.adjust.method = "BH"))
View(ccb$p.value)
class(ccb$p.value)

