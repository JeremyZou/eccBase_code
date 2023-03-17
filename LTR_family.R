#开始组建环境
#rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling","IRanges",
         "GenomicRanges","rtracklayer")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入TE的数据
repeat_infor="/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/mm10_rmsk_TE.gtf"
rpeat<-read.delim(repeat_infor,header=F,stringsAsFactors = F)
colnames(rpeat)<-c("chromosome","a","b","start","end","vaue","strand","oth","feature")
rpeat<-rpeat[which(grepl("LTR", rpeat$feature)),c(1,4,5,9)]
rpeat<-makeGRangesFromDataFrame(rpeat,keep.extra.columns=T)
#开始整理数据：与eccDNA进行findOverlap
f_all <- list.files(path = "./tempTissues_mmu",pattern = ".txt", full.names = T)
d_LTR_family <- data.frame()
for (i in seq(length(list.files(path = "./tempTissues_mmu",pattern = ".txt", full.names = T)))) {
  print(str_split(basename(f_all[i]),pattern = "\\.")[[1]][1])
  temp         <- read.csv(f_all[i], header = F, sep = "\t", stringsAsFactors = F)
  temp  <- temp[,c(1,2)]
  temp  <- na.omit(temp)
  temp  <- separate(data = temp, col = V2, into = c("chromosome", "start","end"), sep = "[-|:]")
  temp[,c(3)] <- as.numeric(temp[,c(3)])
  temp[,c(4)] <- as.numeric(temp[,c(4)])
  temp  <- makeGRangesFromDataFrame(temp, keep.extra.columns = T)
  EDG   <- GenomicRanges::findOverlaps(temp, rpeat, ignore.strand = T)
  data  <- sapply(str_split(rpeat[subjectHits(EDG)]$feature, pattern = ";"),function(x){x[3]})
  data  <- sapply(str_split(data, pattern = " family_id "),function(x){x[2]})
  data  <- table(data)
  data  <- round(prop.table(data),3)
  data  <- as.data.frame(data)
  data$mark <- str_split(basename(f_all[i]),pattern = "\\.")[[1]][1]
  d_LTR_family <- rbind(d_LTR_family, data)
}
#开始作图：与eccDNA进行findOverlap
p_1 <- ggplot(data = d_LTR_family) +
       geom_bar(mapping = aes(x=mark,weight=Freq, fill=data),position = "stack") +
       theme_classic() +
       ylab(label = "Frequency") +
       labs(fill = "LTR repFamily") +
       theme(axis.text.x = element_text(angle = 30, hjust = 0.9),
             plot.margin = unit(c(0.2,0.2,0.2,2),"cm"),
             axis.text = element_text(family = "serif", size = 12, color = "black"),
             axis.title.x = element_blank(),
             axis.title.y = element_text(family = "serif", size = 12, color = "black"),
             legend.text = element_text(family = "serif", size = 12, color = "black"),
             legend.title = element_text(family = "serif", size = 12, color = "black")) +
       scale_fill_nejm()
p_1



#结论：repFamily差不多



