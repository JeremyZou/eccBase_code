#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling","IRanges",
         "GenomicRanges","rtracklayer","data.table")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入TE的数据
repeat_infor = "/Reference/zyw_reference/zyw_reference/table_repeat/Te_transcript/mm10_rmsk_TE.gtf"
rpeat <- read.delim(repeat_infor,header=F,stringsAsFactors = F)
colnames(rpeat) <- c("chromosome","a","b","start","end","vaue","strand","oth","feature")
rpeat$len <- rpeat$end - rpeat$start
rpeat <- rpeat[which(grepl("LTR", rpeat$feature)),c(1,4,5,9,10)]
rpeat <- makeGRangesFromDataFrame(rpeat, keep.extra.columns = T)
#开始整理数据：与eccDNA进行findOverlap
f_all <- list.files(path = "./tempTissues_mmu", pattern = ".txt", full.names = T)
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
for (k in c("ERV1","ERVK","ERVL","ERVL-MaLR","Gypsy","LTR")) {
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
data <- data.frame(mark <- str_split(basename(f_all[i]), pattern = "\\.")[[1]][1],
                   frac <- round(sum(qwe$frac.sum)/length(temp$V4),4),
                   sd   <- round(sd(qwe$frac.sum),4),
                   repFamily <- sapply(k, function(x){if(x==""){x="classLTR"} else if (x==x){x=x}}))
colnames(data) <- c("mark","frac","sd","repFamily")
repFam_frac <- rbind(repFam_frac,data)
}
}
#开始整理数据：计算repClass LTR在各个组织eccDNA的长度占比
aab <- repFam_frac$mark
for (k in unique(aab)) {
  print(k)
  sst  <- repFam_frac[which(repFam_frac$mark==k),]
  data <- data.frame(a<-k,b<-sum(sst$frac),c<-"repClass_LTR")
  colnames(data) <- colnames(sst)
  repFam_frac <- rbind(repFam_frac,data)
}

#开始整理数据：计算LTR各个repFamily的总长度在mm10（chr1-19,X,Y）的占比
for (k in c("ERV1","ERVK","ERVL","ERVL-MaLR","Gypsy","LTR")) {
  print(k)
  inter  <- rpeat[which(grepl(pattern = paste0("family_id ",k,";"), 
                             x = rpeat$feature))]
  inter  <- as.data.frame(inter)
  ran    <- c(paste0("chr",seq(1,19)),"chrX","chrY")
  inter  <- inter[which(inter$seqnames %in% ran),]
  fracAl <- sum(inter$len)/2725521370
  print(fracAl)
}
for (k in c("LTR_repClass")){
inter <- rpeat[which(grepl(pattern = "family_id ", x = rpeat$feature))]
inter <- as.data.frame(inter)
ran   <- c(paste0("chr",seq(1,19)),"chrX","chrY")
inter <- inter[which(inter$seqnames %in% ran),]
fracAl <- sum(inter$len)/2725521370
print(paste0("LTR_repClass  ",fracAl))
}
#LTR repFamily、LTR repClass长度占比结果如下
#[1] "ERV1"
#[1] 0.01141918
#[1] "ERVK"
#[1] 0.04808627
#[1] "ERVL"
#[1] 0.01195779
#[1] "ERVL-MaLR"
#[1] 0.04496413
#[1] "Gypsy"
#[1] 0.0001670968
#[1] "LTR"
#[1] 0.0001240779
#"LTR" repClass
#[1] 0.1167185


#开始作图：repFamily在mmu各个组织的比较
sst <- repFam_frac[which(repFam_frac$repFamily == "repClass_LTR"),]
sst <- sst[order(sst$frac, decreasing = T),]
repF_1 <- repFam_frac[which(repFam_frac$repFamily == "repClass_LTR"),]
repF_2 <- repFam_frac[which(repFam_frac$repFamily != "repClass_LTR"),]
repF_1$mark <- factor(repF_1$mark, levels = sst$mark)
repF_2$mark <- factor(repF_2$mark, levels = sst$mark)
ggplot(data = repF_1) +
  geom_bar(aes(x=repFamily,weight=frac,fill=mark),position = "dodge") +
  scale_fill_d3(palette = "category20") +
  ylab(label = "Fraction") +
  labs(fill="Tissues in mmu") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0,0.12,by=0.03)) +
  guides(fill=guide_legend(ncol=1)) +
  theme(axis.title = element_text(family = "serif", size = 12, color="black"),
        legend.text = element_text(family = "serif", size = 10, color="black"),
        axis.text = element_text(family = "serif", size = 12, color="black"),
        legend.title = element_text(family = "serif", size = 10, color="black"),
        legend.key.size = unit(0.2,"cm"))# +
  annotate("segment",x=6.5,xend = 7.5,y=0.1167185,yend = 0.1167185,
           lty=2,color="red")





p_1 <- ggplot(data = repFam_frac) +
       geom_bar(aes(x=repFamily,weight=frac,fill=mark),position = "dodge") +
       scale_fill_d3(palette = "category20") +
       ylab(label = "Fraction") +
       labs(fill="Tissues in mmu") +
       theme_classic() +
       scale_y_continuous(breaks = seq(0,0.12,by=0.03)) +
       guides(fill=guide_legend(ncol=1)) +
       theme(axis.title = element_text(family = "serif", size = 12, color="black"),
             legend.text = element_text(family = "serif", size = 10, color="black"),
             axis.text = element_text(family = "serif", size = 12, color="black"),
             legend.title = element_text(family = "serif", size = 10, color="black"),
             legend.key.size = unit(0.2,"cm")) +
       annotate("segment",x=0.5,xend = 1.5,y=0.01141918,yend = 0.01141918,
                lty=2,color="red") +
       annotate("segment",x=1.5,xend = 2.5,y=0.04808627,yend = 0.04808627,
                lty=2,color="red") +
       annotate("segment",x=2.5,xend = 3.5,y=0.01195779,yend = 0.01195779,
                lty=2,color="red") +
       annotate("segment",x=3.5,xend = 4.5,y=0.04496413,yend = 0.04496413,
                lty=2,color="red") +
       annotate("segment",x=4.5,xend = 5.4,y=0.0001670968,yend = 0.0001670968,
                lty=2,color="red",size=0.2) +
       annotate("segment",x=5.6,xend = 6.5,y=0.0001240779,yend = 0.0001240779,
                lty=2,color="red",size=0.2) +
       annotate("segment",x=6.5,xend = 7.5,y=0.1167185,yend = 0.1167185,
                lty=2,color="red") +
       annotate("text",label="Red dash line mean \naverage level of LTR in mm10",
                x=3.5,y=0.1,color="red", family="serif",hjust=0) 
p_1
ggsave(p_1, filename = "./tempTissues_mmu/LTR_famFrac.jpg", dpi = 800, 
       units = "cm", width = 24, height=10)






