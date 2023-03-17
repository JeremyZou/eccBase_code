#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling","IRanges",
         "GenomicRanges","rtracklayer","data.table","patchwork","tidyr")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入TE的数据
repeat_infor = "/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/mm10_rmsk_TE.gtf"
rpeat <- read.delim(repeat_infor,header=F,stringsAsFactors = F)
colnames(rpeat) <- c("chromosome","a","b","start","end","vaue","strand","oth","feature")
rpeat$len <- rpeat$end - rpeat$start
rpeat <- rpeat[which(grepl("LTR", rpeat$feature)), c(1, 4, 5, 9, 10)]
rpeat <- makeGRangesFromDataFrame(rpeat, keep.extra.columns = T)
#开始整理数据：与eccDNA进行findOverlap
f_all <- list.files(path = "./tempTissues_mmu", pattern = ".txt", full.names = T)
aa = ".*"#画p_1
aa=c("ERV1","ERVK","ERVL","ERVL-MaLR","Gypsy","LTR")#画p_2

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
    qwe$frac <- round(qwe$ltrL/qwe$eccL, 3)
    qwe <- as.data.table(qwe)
    qwe <- qwe[,.(frac.sum = sum(frac)),by=eccN]
    data <- data.frame(mark <- rep(str_split(basename(f_all[i]), pattern = "\\.")[[1]][1]),
                       frac <- round(sum(qwe$frac.sum)/length(temp$V4),4),
                       sd   <- round(sd(c(qwe$frac.sum,rep(0,length(temp$V4)-length(qwe$frac.sum))))/sqrt(length(temp$V4)),4),
                       repFamily <- k)
    colnames(data) <- c("mark","frac","sd","repFamily")
    repFam_frac <- rbind(repFam_frac,data)
  }
}


#开始作图：repClass LTR在小鼠组织eccDNA的平均占比
repCla_frac <- repFam_frac
repCla_frac$repFamily <- "LTR"
repCla_frac <- repCla_frac[order(repCla_frac$frac, decreasing = T),]
repCla_frac$mark <- factor(repCla_frac$mark, levels = unique(repCla_frac$mark))
p_1 <- ggplot(data = repCla_frac,aes(x=repFamily,weight=frac,fill=mark)) +
       geom_bar(position = "dodge") +
       geom_errorbar(aes(ymin=frac, ymax=frac+sd), position = position_dodge(0.9), width=0.5, size=0.2) +
       scale_fill_d3(palette = "category20") +
       ylab(label = "Fraction") +
       xlab(label = "repClass") +
       labs(fill="Tissues in mmu") +
       theme_classic() +
       scale_x_discrete(expand = c(0.5, 0)) +
       scale_y_continuous(expand = c(0.02, 0),limits = c(0, 0.13), breaks = seq(0, 0.15, by = 0.03)) +
       #guides(fill = guide_legend(ncol = 1)) +
       guides(fill="none") +
       theme(axis.title = element_text(family = "serif", size = 12, color="black"),
             legend.text = element_text(family = "serif", size = 10, color="black"),
             axis.text = element_text(family = "serif", size = 10, color="black"),
             legend.title = element_text(family = "serif", size = 10, color="black"),
             legend.key.size = unit(0.2,"cm")) +
       annotate("segment",x=0.55,xend = 1.45,y=0.1167185,yend = 0.1167185,
              lty=2,color="red")
p_1


#开始作图：repFamily LTR在小鼠组织eccDNA的平均占比
repFam_frac$mark <- factor(repFam_frac$mark, levels = repCla_frac$mark)
p_2 <- ggplot(data = repFam_frac,aes(x=repFamily,weight=frac,fill=mark)) +
  geom_bar(position = "dodge") +
  geom_errorbar(aes(ymin=frac, ymax=frac+sd), position = position_dodge(0.9), width=0.5, size=0.2) +
  scale_fill_d3(palette = "category20") +
  ylab(label = "Length fraction in eccDNA") +
  labs(fill="Tissues in mmu") +
  theme_classic() +
  scale_x_discrete(expand = c(0.1, 0)) +
  scale_y_continuous(expand = c(0.02, 0),limits = c(0, 0.13), breaks = seq(0, 0.15, by = 0.03)) +
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
        annotate("segment",x=0.55,xend = 1.45,y=0.01141918,yend = 0.01141918,
                 lty=2,color="red") +
        annotate("segment",x=1.55,xend = 2.45,y=0.04808627,yend = 0.04808627,
                 lty=2,color="red") +
        annotate("segment",x=2.55,xend = 3.45,y=0.01195779,yend = 0.01195779,
                 lty=2,color="red") +
        annotate("segment",x=3.55,xend = 4.45,y=0.04496413,yend = 0.04496413,
                 lty=2,color="red") +
        annotate("segment",x=4.55,xend = 5.45,y=0.0001670968,yend = 0.0001670968,
                 lty=2,color="red",size=0.2) +
        annotate("segment",x=5.55,xend = 6.45,y=0.0001240779,yend = 0.0001240779,
                 lty=2,color="red",size=0.2)
p_2

myplot_ltr <- p_1 + p_2 + plot_layout(ncol = 2, widths = c(1,7))
myplot_ltr
getwd()
#ggsave(myplot_ltr,filename = "tempTissues_mmu/1_ltr_repClassrepFamily.jpg",dpi = 800, units = "cm", width = 24, height = 7)
ggsave(myplot_ltr,filename = "tempTissues_mmu/1_ltr_repClassrepFamily.pdf",units = "cm", width = 20, height = 7)



