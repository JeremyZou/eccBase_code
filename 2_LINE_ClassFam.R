#开始组建环境
#rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling","IRanges",
         "GenomicRanges","rtracklayer","data.table","patchwork")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入TE的数据
repeat_infor = "/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/mm10_rmsk_TE.gtf"
rpeat <- read.delim(repeat_infor, header=F, stringsAsFactors = F)
colnames(rpeat) <- c("chromosome","a","b","start","end","vaue","strand","oth","feature")
rpeat$len <- rpeat$end - rpeat$start
rpeat <- rpeat[which(grepl("LINE", rpeat$feature)), c(1,4,5,9,10)]
rpeat <- makeGRangesFromDataFrame(rpeat, keep.extra.columns = T)
#开始整理数据：与eccDNA进行findOverlap
f_all <- list.files(path = "./tempTissues_mmu", pattern = ".txt", full.names = T)
repFam_fracLINE <- data.frame()
aa = ".*"
aa = c("CR1","Dong-R4","L1","L2","Penelope","RTE-BovB","RTE-X")
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
                                      rpeat[which(grepl(pattern = paste0("family_id ",k,";"), x = rpeat$feature))], 
                                      ignore.strand = T)
    inter <- as.data.frame(inter)
    inter$feature <- inter$width - 1
    inter <- makeGRangesFromDataFrame(inter, keep.extra.columns = T)
    EDG   <- GenomicRanges::findOverlaps(temp, inter, ignore.strand = T)
    qwe   <- data.frame(eccN <- temp[queryHits(EDG)]$V1,
                        eccL <- temp[queryHits(EDG)]$V4,
                        ltrL <- inter[subjectHits(EDG)]$feature)
    colnames(qwe) <- c("eccN","eccL","ltrL")
    qwe$frac <- round(qwe$ltrL/qwe$eccL, 4)
    qwe <- as.data.table(qwe)
    qwe <- qwe[,.(frac.sum = sum(frac)),by=eccN]
    data <- data.frame(mark <- rep(str_split(basename(f_all[i]), pattern = "\\.")[[1]][1]),
                       frac <- round(sum(qwe$frac.sum)/length(temp$V4),4),
                       sd   <- round(sd(c(qwe$frac.sum,rep(0,length(temp$V4)-length(qwe$frac.sum))))/sqrt(length(temp$V4)),4),
                       repFamily <- k)
    colnames(data) <- c("mark","frac","sd","repFamily")
    repFam_fracLINE <- rbind(repFam_fracLINE,data)
  }
}

#开始整理数据：计算LINE各个repFamily的总长度在mm10（chr1-19,X,Y）的占比
for (k in c("CR1","Dong-R4","L1","L2","Penelope","RTE-BovB","RTE-X")) {
  print(k)
  inter <- rpeat[which(grepl(pattern = paste0("family_id ",k,";"), 
                             x = rpeat$feature))]
  inter <- as.data.frame(inter)
  ran   <- c(paste0("chr",seq(1,19)),"chrX","chrY")
  inter <- inter[which(inter$seqnames %in% ran),]
  fracAl <- sum(inter$len)/2725521370
  print(fracAl)
}
inter <- rpeat[which(grepl(pattern = "family_id ", x = rpeat$feature))]
inter <- as.data.frame(inter)
ran   <- c(paste0("chr",seq(1,19)),"chrX","chrY")
inter <- inter[which(inter$seqnames %in% ran),]
fracAl <- sum(inter$len)/2725521370
print(fracAl)
#长度占比结果如下
#####,,,~~~~"CR1"
#####,,,~~~~0.0006536364
#####,,,~~~~"Dong-R4"
#####,,,~~~~0.000008935171
#####,,,~~~~"L1"
#####,,,~~~~0.195713
#####,,,~~~~"L2"
#####,,,~~~~0.003943057
#####,,,~~~~"Penelope"
#####,,,~~~~0.000002981448
#####,,,~~~~"RTE-BovB"
#####,,,~~~~0.000007968751
#####,,,~~~~"RTE-X"
#####,,,~~~~0.0001105447
#####,,,~~~~repClass LINE
#####,,,~~~~0.2004401



#开始作图：repClass LINE在小鼠组织eccDNA的平均占比
repCla_fracLINE <- repFam_fracLINE
repCla_fracLINE$repFamily <- "LINE"
repCla_fracLINE <- repCla_fracLINE[order(repCla_fracLINE$frac, decreasing = T),]
repCla_fracLINE$mark <- factor(repCla_fracLINE$mark, levels = repCla_fracLINE$mark)
p1_LINE <- ggplot(data = repCla_fracLINE,aes(x=repFamily,weight=frac,fill=mark)) +
           geom_bar(position = "dodge") +
           geom_errorbar(aes(ymin=frac, ymax=frac+sd), position = position_dodge(0.9), width=0.5, size=0.2) +
           scale_fill_d3(palette = "category20") +
           ylab(label = "Fraction") +
           xlab(label = "repClass") +
           labs(fill="Tissues in mmu") +
           theme_classic() +
           scale_x_discrete(expand = c(0.5, 0)) +
           scale_y_continuous(expand = c(0.02, 0),limits = c(0, 0.22), breaks = seq(0, 0.22, by = 0.05)) +
           #guides(fill = guide_legend(ncol = 1)) +
           guides(fill="none") +
           theme(axis.title = element_text(family = "serif", size = 12, color="black"),
                 legend.text = element_text(family = "serif", size = 10, color="black"),
                 axis.text = element_text(family = "serif", size = 10, color="black"),
                 legend.title = element_text(family = "serif", size = 10, color="black"),
                 legend.key.size = unit(0.2,"cm")) +
           annotate("segment",x=0.55,xend = 1.45,y=0.2004401,yend = 0.2004401, lty=2,color="red")
p1_LINE



#开始作图：repFamily LINE在小鼠组织eccDNA的平均占比
repFam_fracLINE$mark <- factor(repFam_fracLINE$mark, levels = repCla_fracLINE$mark)
p2_LINE <- ggplot(data = repFam_fracLINE,aes(x = repFamily, weight = frac, fill = mark)) +
  geom_bar(position = "dodge") +
  geom_errorbar(aes(ymin=frac, ymax=frac+sd), position = position_dodge(0.9), width=0.5, size=0.2) +
  scale_fill_d3(palette = "category20") +
  ylab(label = "Length fraction in eccDNA") +
  labs(fill="Tissues in mmu") +
  theme_classic() +
  scale_x_discrete(expand = c(0.1, 0)) +
  scale_y_continuous(expand = c(0.02, 0),limits = c(0, 0.22), breaks = seq(0, 0.22, by = 0.05)) +
  guides(fill = guide_legend(ncol = 1)) +
  #guides(fill="none") +
  theme(axis.title = element_text(family = "serif", size = 12, color="black"),
        legend.text = element_text(family = "serif", size = 10, color="black"),
        axis.text.x = element_text(family = "serif", size = 10, color="black"),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        #plot.margin = unit(rep(0.2,4),"cm"),
        axis.text.y = element_blank(),
        legend.title = element_text(family = "serif", size = 12, color="black"),
        legend.key.size = unit(0.2,"cm")) +
  annotate("segment",x=0.55,xend = 1.45,y=0.0006536364,yend = 0.0006536364,
           lty=2,color="red",size=0.1) +
  annotate("segment",x=1.55,xend = 2.45,y=0.000008935171,yend = 0.000008935171,
           lty=2,color="red",size=0.1) +
  annotate("segment",x=2.55,xend = 3.45,y=0.195713,yend = 0.195713,
           lty=2,color="red") +
  annotate("segment",x=3.55,xend = 4.45,y=0.003943057,yend = 0.003943057,
           lty=2,color="red") +
  annotate("segment",x=4.55,xend = 5.45,y=0.000002981448,yend = 0.000002981448,
           lty=2,color="red",size=0.1) +
  annotate("segment",x=5.55,xend = 6.45,y=0.000007968751,yend = 0.000007968751,
           lty=2,color="red",size=0.1) +
  annotate("segment",x=6.55,xend = 7.45,y=0.0001105447,yend = 0.0001105447,
           lty=2,color="red",size=0.1)
p2_LINE

myplot_LINE <- p1_LINE + p2_LINE + plot_layout(ncol = 2, widths = c(1,7))
myplot_LINE
getwd()
#ggsave(myplot_LINE,filename = "tempTissues_mmu/2_LINE_repClassrepFamily.jpg",dpi = 800, units = "cm", width = 20, height = 7)
ggsave(myplot_LINE,filename = "tempTissues_mmu/2_LINE_repClassrepFamily.pdf",units = "cm", width = 20, height = 7)














#开始整理数据：计算repClass LTR在各个组织eccDNA的长度占比
#aab <- repFam_fracLINE$mark
#for (k in unique(aab)) {
#  print(k)
#  sst  <- repFam_fracLINE[which(repFam_fracLINE$mark==k),]
#  data <- data.frame(a<-k, b<-sum(sst$frac), c<-"repClass_LINE")
#  colnames(data) <- colnames(sst)
#  repFam_fracLINE <- rbind(repFam_fracLINE,data)
#}


