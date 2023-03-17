#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap")
for (i in pks){library(i, character.only = T)}
rm(i,pks)


#开始整理数据：智人的相关数据
ha <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/hsa2.txt", sep="\t", stringsAsFactors = F, header = T)
hb <- cbind(ha[, grepl("count",colnames(ha))], ha[, c(1, 5, 6, 48)])
uniqA <- as.data.frame(apply((as.data.frame(table(hb$repeat_count))),2,function(x){
  as.integer(x)
}))
uniqA$Lo <-log2(uniqA$Freq +1)
#开始整理数据：小鼠的相关数据
ma <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/mmu.txt",
               sep="\t",stringsAsFactors = F,header = T)
mb <- cbind(ma[,grepl("count",colnames(ma))],ma[,c(1,5,6,48)])
uniqM <- as.data.frame(apply((as.data.frame(table(mb$repeat_count))),2,function(x){
  as.integer(x)
}))
uniqM$Lo <-log2(uniqM$Freq +1)
#开始整理数据：智人的Healthy、Cancer组别数据
he <- hb %>% dplyr::filter(Disease == "Healthy")
ce <- hb %>% dplyr::filter(Disease != "Healthy")
uniqH <- as.data.frame(apply((as.data.frame(table(he$repeat_count))),2,function(x){
  as.integer(x)
}))
uniqH$Lo <-log2(uniqH$Freq +1)
uniqC <- as.data.frame(apply((as.data.frame(table(ce$repeat_count))),2,function(x){
  as.integer(x)
}))
uniqC$Lo <-log2(uniqC$Freq +1)


#开始作图：智人的柱状图
p_1 <- ggplot(data = uniqA) +
       theme_classic() +
       theme(axis.title = element_text(family = "serif", size = 12),
             axis.text  = element_text(family = "serif", size = 12)) +
       geom_bar(aes(x=Var1,y=Lo), stat = "identity",width = 0.2,color="indianred3",
                fill="indianred3") +
       labs(y="Log2(Number of eccNDAs + 1)",x="Number of eccDNA repeats") +
       scale_y_continuous(limits = c(0,21),breaks = seq(0,20, by=2),
                          expand = c(0.015,0)) +
       scale_x_continuous(breaks = seq(0,420000, by=100000),expand = c(0.015,0))
p_1
ggsave(p_1, filename = "repeat_hsa.jpg", device = "jpeg", dpi = 800)
#开始作图：小鼠的柱状图
p_2 <- ggplot(data = uniqM) +
       theme_classic() +
       theme(axis.title = element_text(family = "serif",size = 12),
             axis.text  = element_text(family = "serif",size = 12)) +
       geom_bar(aes(x=Var1,y=Lo), stat = "identity",width = 0.2,color="indianred3",
                fill="indianred3")+
       scale_y_continuous(limits = c(-0.2,20),breaks = seq(0,20,by=4),expand = c(0.015,0)) +
       scale_x_continuous(breaks = seq(0,16,by=2),expand = c(0.015,0)) +
       labs(y="Log2(Number of eccNDAs + 1)",x="Number of eccDNA repeats") 
p_2
ggsave(p_2, filename = "repeat_mmu.jpg", device = "jpeg", dpi = 800)
#开始作图：智人的柱状图 + 频率折线图
jpeg(filename = "repeat_histogram_line_hsa.jpg", res = 800,units = "cm", width = 8, height = 4)
p_3 <- ggplot() +
       theme_classic() +
       geom_bar(data = uniqA, aes(x=Var1,y=Lo), stat = "identity",width = 0.2,
                color="indianred3",alpha=1) +
       scale_x_continuous(expand = c(0.02, 0), breaks = c(seq(0,420000,by=100000))) +
       scale_y_continuous(limits = c(0,20),expand = c(0.035, 0), breaks = seq(0,20,by=4),
                          sec.axis = sec_axis(trans = ~.,breaks = seq(0,20,by=4),
                                              labels = c(0,paste0(seq(0,20,by=4)[2:6]*3,"%")),
                                              name = "Frequency")) +
       geom_freqpoly(data = hb, aes(repeat_count, y = ..count../736775*100/3), 
                     binwidth = 1, color="#1F77B4FF", size=1, alpha=1) +
       theme(axis.title.x       = element_text(family = "serif",color = "black",size = 10),
             axis.title.y.left  = element_text(family = "serif",color = "indianred3",size = 10),
             axis.title.y.right = element_text(family = "serif",color = "#1F77B4FF",size  = 10,
                                               angle = 90,vjust = 1),
             axis.text.y.left   = element_text(size = 10,family = "serif",color = "indianred3"),
             axis.text.y.right  = element_text(size = 10,family = "serif",color = "#1F77B4FF"),
             axis.line.y.left   = element_line(colour = "indianred3",size = 0.8),
             axis.line.y.right  = element_line(colour = "#1F77B4FF",size = 0.8),
             axis.ticks.y.left  = element_line(colour = "indianred3",size = 0.8),
             axis.ticks.y.right = element_line(colour = "#1F77B4FF",size = 0.8)) +
       xlab(expression(Number~of~eccDNA~repeats)) +ylab("Log2(Number of\n eccDNAs + 1)")
p_3
ggsave(plot = p_3, filename = "repeat_histogram_line_hsa.pdf", units = "cm",width = 8, height = 4)
dev.off()
#开始作图：小鼠的柱状图 + 频率折线图
jpeg(filename = "repeat_histogram_line_mmu.jpg", res = 800,units = "cm", width = 8, height = 4)
p_4 <- ggplot() +
       theme_classic() +
       geom_bar(data = uniqM, aes(x=Var1,y=Lo), stat = "identity",width = 0.2,
                color="indianred3",fill="indianred3",alpha=1) +
       scale_x_continuous(expand = c(-0.058, 0), breaks = c(seq(0,16,by=4))) +
       scale_y_continuous(limits = c(0,20),expand = c(0.02, 0), breaks = seq(0,20,by=4),
                          sec.axis = sec_axis(trans = ~.,breaks = seq(0,20,by=4),
                                              labels = c(0,paste0(seq(0,20,by=4)[2:6]*3,"%")),
                                              name = "Frequency")) +
       geom_freqpoly(data = mb, aes(repeat_count, y = ..count../481831*100/3), 
                     binwidth = 1, color="#1F77B4FF", size=1, alpha=1) +
       theme(axis.title.x       = element_text(family = "serif",color = "black",size = 10),
             axis.title.y.left  = element_text(family = "serif",color = "indianred3",size = 10),
             axis.title.y.right = element_text(family = "serif",color = "#1F77B4FF",size  = 10,
                                               angle = 90,vjust = 1),
             axis.text.y.left   = element_text(size = 10,family = "serif",color = "indianred3"),
             axis.text.y.right  = element_text(size = 10,family = "serif",color = "#1F77B4FF"),
             axis.line.y.left   = element_line(colour = "indianred3",size = 0.8),
             axis.line.y.right  = element_line(colour = "#1F77B4FF", size = 0.8),
             axis.ticks.y.left  = element_line(colour = "indianred3",size = 0.8),
             axis.ticks.y.right = element_line(colour = "#1F77B4FF", size = 0.8)) +
       xlab(expression(Number~of~eccDNA~repeats)) +ylab("Log2(Number of\n eccDNAs + 1)")
p_4
ggsave(plot = p_4, filename = "repeat_histogram_line_mmu.pdf",units = "cm", width = 8, height = 4)
dev.off()


#开始统计数据：关于智人的所有"frac"字段的统计（不为0的eccDNA数量占比）
for (i in seq(length(colnames(hb))-3)) {
  print(colnames(hb[i]))
  tempC <- length(hb[which(hb[,i]>0),][,i])/length(hb[,1])
  print(tempC)
}




#其他相关代码
#da[,length(colnames(da))+1] <- apply(da[,c(1,2)], 1, function(x) {
#  #inds <- match(x[1], unicq$Var1)
#  inds <- match(as.character(x[1]), as.character(unicq$Var1))
#  ifelse(is.na(inds), x[1], unicq$Freq[inds]) 
#})