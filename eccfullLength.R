#hsa的p_3, mmu的p_4是最终图
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_length/")
ha <- read.csv(file = "../raw/hsa2.txt",sep = "\t",stringsAsFactors = F,header = T)
#开始整理数据
aa_2h <- ha[,c("Length","eccDNA_ID")]
#1、for using geom_bar() to plot
dataBar_1 <- data.frame(ca1<-length(aa_2h[which(aa_2h$Length < 100),1]), 
                   ca2<-length(aa_2h[which(aa_2h$Length >= 100 & aa_2h$Length < 300 ),1]),
                   ca3<-length(aa_2h[which(aa_2h$Length >= 300 & aa_2h$Length < 600 ),1]),
                   ca4<-length(aa_2h[which(aa_2h$Length >= 600 & aa_2h$Length < 900 ),1]),
                   ca5<-length(aa_2h[which(aa_2h$Length >= 900 & aa_2h$Length < 1200 ),1]),
                   ca6<-length(aa_2h[which(aa_2h$Length >= 1200 & aa_2h$Length < 10000 ),1]),
                   ca7<-length(aa_2h[which(aa_2h$Length >= 10000 & aa_2h$Length < 1000000 ),1]),
                   ca8<-length(aa_2h[which(aa_2h$Length >= 1000000 ),1]))
dataBar_1 <- data.frame(t(dataBar_1))
rownames(dataBar_1) <- c("ca1","ca2","ca3","ca4","ca5","ca6","ca7","ca8")
dataBar_1$bn2       <- rownames(dataBar_1)
#2、for using geom_bar() to plot
dataBar_2    <- as.data.table(aa_2h)
dataBar_2[which(dataBar_2[,1]>=10000 & dataBar_2[,1]<1000000),][,1] <- 12000
dataBar_2[which(dataBar_2[,1]>=1000000),][,1] <- 14000
dataBar_2[which(dataBar_2[,1]>=5000000),][,1] <- 5000000

dataBar_2    <- dataBar_2[, .N, by = c("Length")]
setkey(dataBar_2, N)
#3、for using geom_histogram() to plot
dataHis_1 <- aa_2h
dataHis_1[which(dataHis_1[,1]>=1000000),][,1] <- 1200000

#1、开始作图:bar：柱状图
jpeg(filename = "fullLength_bar_hsa_1.jpg",res = 800,units = "cm",width = 9,height = 6)
p_1 <- ggplot(data = temp) +
       theme_classic() +
       ggtitle(label = "Full length of hsa eccDNA") +
       theme(axis.title = element_text(family = "serif",size = 16),
             axis.text = element_text(family = "serif",size = 16) ,
             plot.title = element_text(family = "serif",size = 16,hjust = 0.5)) +
       geom_bar(aes(y=t.temp.,x=bn2), stat = "identity",width = 0.5,color="indianred3") +
       labs(y="Number of eccNDAs", x="Length of eccDNA") #+
       scale_y_continuous(limits = c(-0.2,12),breaks = seq(0,12,by=1)) +
       scale_x_continuous(limits = c(-0.2,5000000),breaks = c(seq(0,5000000,by=300000)),
                          expand = c(0.015,0))
p_1
dev.off()

#2、开始作图:bar：柱状图
jpeg(filename = "fullLength_bar_hsa_2.jpg", res = 800, units = "cm", width = 9, height = 6)
p_2 <- ggplot(data = dataBar_2) +
       theme_classic() +
       ggtitle(label = "Full length of hsa eccDNA") +
       theme(axis.title = element_text(family = "serif",size = 16),
             axis.text  = element_text(family = "serif",size = 16) ,
             plot.title = element_text(family = "serif",size = 16,hjust = 0.5)) +
       geom_bar(aes(y=N,x=Length), stat = "identity",width = 0.5,color="indianred3") +
       labs(y="Number of eccNDAs", x="Length of eccDNA") #+
       scale_y_continuous(limits = c(-0.2,12),breaks = seq(0,12,by=1)) +
       scale_x_continuous(limits = c(-0.2,5000000),breaks = c(seq(0,5000000,by=300000)),
                          expand = c(0.015,0))
p_2
dev.off()

#3、开始作图:histogram：直方图
jpeg(filename = "fullLength_histogram_hsa.jpg",res = 800,units = "cm",width = 9,height = 6)
p_3 <- ggplot(dataHis_1,aes(x=Length)) +
       geom_histogram(aes(y=log2(..count.. + 1)),binwidth = 1000000,fill="indianred3",alpha=1.2) +
       #geom_freqpoly(aes(Length, y=log2(..count.. + 1)), binwidth = 1000000, color="black", size=0.5) +
       theme_classic() +
       scale_y_continuous(expand = c(0.01,1),breaks = seq(0,20,by=2)) +
       scale_x_continuous(expand = c(0.02,0.01),breaks = seq(0, 2400000000, by = 30000000),
                          labels = seq(0, 2400, by = 30)) +
       theme(axis.title  = element_text(family = "serif",color = "black",size = 12),
             axis.text   = element_text(size = 12,family = "serif",color = "black"),
             axis.line   = element_line(colour = "black",size = 0.8)) +
       xlab("eccDNA size (× 10^6 bp)") +ylab("Log2(Number of eccDNAs + 1)")
p_3
dev.off()

#3、开始作图:histogram：直方图 + 频率折线图
jpeg(filename = "fullLength_histogram_bar_hsa.jpg",res = 800,units = "cm",width = 8,height = 4)
p_3x <- ggplot(dataHis_1,aes(x=Length)) +
  geom_histogram(aes(y=log2(..count.. + 1)),binwidth = 1000000,fill="indianred3",alpha=1.2) +
  #geom_freqpoly(aes(Length, y=log2(..count.. + 1)), binwidth = 1000000, color="black", size=0.5) +
  geom_freqpoly(aes(Length, y=..count../736775*20), binwidth = 1000000, color="#1F77B4FF", size=0.5, alpha=1) +
  theme_classic() +
  scale_y_continuous(expand = c(0.01,0.5),breaks = seq(0,20,by=4),
                     sec.axis = sec_axis(trans = ~.,breaks = seq(0,20,by=4),
                                         labels = c(0,paste0(seq(0,20,by=4)[2:6]*5,"%")),
                                         name = "Frequency")) +
  scale_x_continuous(expand = c(0.02,0.01),breaks = seq(0, 2400000000, by = 30000000),
                     labels = seq(0, 2400, by = 30)) +
  theme(axis.title.x       = element_text(family = "serif",color = "black",size = 10),
        axis.title.y.left  = element_text(family = "serif",color = "indianred3",size = 10),
        axis.title.y.right = element_text(family = "serif",color = "#1F77B4FF",size  = 10,angle = 90,vjust = 2.5),
        axis.text.y.left   = element_text(size = 10,family = "serif",color = "indianred3"),
        axis.text.y.right  = element_text(size = 10,family = "serif",color = "#1F77B4FF"),
        axis.line.y.left   = element_line(colour = "indianred3",size = 0.8),
        axis.line.y.right  = element_line(colour = "#1F77B4FF",size = 0.8),
        axis.ticks.y.left  = element_line(colour = "indianred3",size = 0.8),
        axis.ticks.y.right = element_line(colour = "#1F77B4FF",size = 0.8)) +
  xlab(expression(eccDNA~size~("×"~10^6~bp))) +ylab("Log2(Number of\n eccDNAs + 1)")
p_3x
ggsave(plot = p_3x, filename = "fullLength_histogram_bar_hsa.pdf",units="cm",width = 8,height = 4)
dev.off()
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
#################################################################################################################
ma <- read.csv(file = "../raw/mmu.txt",sep = "\t",stringsAsFactors = F,header = T)
#1、for using geom_histogram() to plot
aa_2m <- ma[,c("Length","eccDNA_ID")]
max(aa_2m[,"Length"])
median(aa_2m[,"Length"])
as.data.frame(table(aa_2m[,"Length"]))[159,]





#2、for using geom_bar() to plot
aa_2mBar <- as.data.table(aa_2m)
aa_2mBar <- aa_2mBar[, .N, by=c("Length")]
setkey(aa_2mBar,N)

#1、开始作图:histogram：直方图
jpeg(filename = "fullLength_histogram_mmu.jpg",res = 800,units = "cm",width = 8,height = 4)
p_4 <- ggplot(aa_2m,aes(x=Length)) +
       geom_histogram(aes(y=log2(..count.. + 1)),binwidth = 1,fill="indianred3") +
       geom_density(aes(Length, y=log2(481831*..density.. + 1)), color="#5050FFFF", size=0.5) +
       theme_classic() +
       scale_y_continuous(expand = c(0.02,0.01),breaks = seq(0,12,by=1),limits = c(0,12)) +
       scale_x_continuous(expand = c(0.02,0.01),breaks = seq(0, 5200, by = 500)) +
       theme(axis.title = element_text(family = "serif",color = "black",size = 12)) +
       theme(axis.text = element_text(size = 12,family = "serif",color = "black")) +
       xlab("eccDNA size (bp)") +ylab("Log2(Number of eccDNAs + 1)")
       #+ coord_cartesian(xlim = c(5000,5200))
p_4
dev.off()
#2、开始作图:histogram：直方图 + 频率折线图
jpeg(filename = "fullLength_histogram_bar_mmu.jpg",res = 800,units = "cm",width = 8,height = 4)
p_4x <- ggplot(aa_2m,aes(x=Length)) +
        geom_histogram(aes(y=log2(..count.. + 1)),binwidth = 1,fill="indianred3",alpha=1) +
        #geom_freqpoly(aes(Length, y=log2(..count.. + 1)), binwidth = 1000000, color="black", size=0.5) +
        geom_freqpoly(aes(Length, y=..count../481831*100*10), binwidth = 1, color="#1F77B4FF", size=0.5, alpha=1) +
        theme_classic() +
        scale_y_continuous(expand = c(0.01,0.4),breaks = seq(0,12,by=2),
                           sec.axis = sec_axis(trans = ~.,breaks = seq(0,12,by=2),
                                               labels = c(0,paste0(seq(0.2,1,by=0.2)/10*100/10,"%")," "),
                                               name = "Frequency")) +
        scale_x_continuous(expand = c(0.02,0.01),breaks = seq(0, 5200, by = 1000)) +
        theme(axis.title.x       = element_text(family = "serif",color = "black",size = 10),
              axis.title.y.left  = element_text(family = "serif",color = "indianred3",size  = 10),
              axis.title.y.right = element_text(family = "serif",color = "#1F77B4FF", size  = 10,angle = 90,vjust = 0.5),
              axis.text.y.left   = element_text(size = 10,family = "serif",color = "indianred3"),
              axis.text.y.right  = element_text(size = 10,family = "serif",color = "#1F77B4FF"),
              axis.line.y.left   = element_line(colour = "indianred3",size = 0.8),
              axis.line.y.right  = element_line(colour = "#1F77B4FF",size = 0.8),
              axis.ticks.y.left  = element_line(colour = "indianred3",size = 0.8),
              axis.ticks.y.right = element_line(colour = "#1F77B4FF",size = 0.8)) +
        #coord_cartesian(xlim = c(4800,5200)) +
        xlab(expression(eccDNA~size~(bp))) +ylab("Log2(Number of\n eccDNAs + 1)")
p_4x
getwd()
ggsave(plot = p_4x, filename = "fullLength_histogram_bar_mmu.pdf",units="cm",width = 8,height = 4)
dev.off()
#3、开始作图:bar：柱状图
jpeg(filename = "fullLength_bar_hsa_1.jpg",res = 800,units = "cm",width = 9,height = 6)
p_5 <- ggplot(data = aa_2mBar) +
       theme_classic() +
       ggtitle(label = "Full length of hsa eccDNA") +
       theme(axis.title = element_text(family = "serif", size = 16),
             axis.text  = element_text(family = "serif", size = 16) ,
             plot.title = element_text(family = "serif", size = 16, hjust = 0.5)) +
       geom_bar(aes(y=N,x=Length), stat = "identity", width = 0.5,color="indianred3") +
       labs(y="Number of eccNDAs", x="Length of eccDNA") #+
       scale_y_continuous(limits = c(-0.2,12), breaks = seq(0,12,by=1)) +
       scale_x_continuous(limits = c(-0.2,5000000),breaks = c(seq(0,5000000,by=300000)),
                          expand = c(0.015,0))
p_5
dev.off()

