#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_length")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap")
for (i in pks){library(i, character.only = T)}
rm(i,pks)


#开始整理数据：智人的Healthy、Cancer组别的Length数据，并限定范围
ha   <- read.csv(file = "../raw/hsa2.txt", sep = "\t", stringsAsFactors = F, header = T)
length(ha[which(ha$Length<51000000 & ha$Length>50000000),"Length"])
max(ha[,"Length"])

ha_h <- ha[grep(pattern = "Healthy", ha[,"Disease"]), c("Disease", "eccDNA_ID", "Length")]
ha_h <- ha_h[which(ha_h$Length < 3000),]
ha_c <- ha[grep(pattern = "[^Healthy]", ha[,"Disease"]), c("Disease", "eccDNA_ID", "Length")]
ha_c <- ha_c[which(ha_c$Length < 3000),]
#开始整理数据：智人的不同检测方法组别的Length数据，并限定范围
ha_cm <- ha[grep(pattern = "[^Healthy]", x = ha[,"Disease"]), c("Disease", "eccDNA_ID", "Length", "Identify_method")]
ha_m1 <- ha_cm[grep(pattern = "Island",  x = ha_cm[, "Identify_method"]), c("Disease", "eccDNA_ID", "Length")]
ha_m1 <- na.omit(ha_m1); ha_m1 <- ha_m1[which(ha_m1$Length <3000),]
ha_m2 <- ha_cm[grep(pattern = "Predict|Manually", x = ha_cm[, "Identify_method"]), c("Disease", "eccDNA_ID", "Length")]
ha_m2 <- na.omit(ha_m2); ha_m2 <- ha_m2[which(ha_m2$Length <3000),]
#开始整理数据：小鼠的Length数据，并限定范围
ma    <- read.csv(file = "../raw/mmu.txt", sep = "\t", stringsAsFactors = F, header = T)
ma_dA <- ma[grep(pattern = "Healthy", x = ma[, "Disease"]), c("Disease", "eccDNA_ID", "Length")]
ma_dA <- ma_dA[which(ma_dA$Length < 3000),]


#开始作图：智人的Healthy组别的直方图
jpeg(filename = "Healthy_hsa_red.jpg", res = 800, units = "cm", width = 9, height = 6)
p_1 <- ggplot(ha_h) + 
       geom_histogram(aes(x = Length), binwidth = 1, fill = "royalblue4") +
       theme_classic() +
       #theme_get() +
       scale_y_continuous(limits = c(0, 1300), expand = c(0.03, 0.01), breaks = seq(0, 1300, by = 400)) +
       scale_x_continuous(expand = c(0.04, 0.01), breaks = seq(0, 3100, by = 500)) +
       theme(axis.title = element_text(family = "serif", color = "black", size = 12),
             axis.text  = element_text(family = "serif", color = "black", size = 12),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "white"),
             panel.grid.minor = element_line(colour = "white")) +
       xlab("eccDNA size (bp)") + ylab("Number of eccDNA")
p_1
getwd()
dir()
ggsave(plot = p_1, filename = "Healthy_hsa_under3000bp.pdf", units = "cm", width = 18, height = 12)
dev.off()
#开始作图：智人的Cancer组别的直方图
jpeg(filename = "Cancer_hsa_red.jpg", res = 800, units = "cm", width = 9, height = 6.1)
p_2 <- ggplot(ha_c) +
       geom_histogram(aes(x = Length), binwidth = 1, fill = "royalblue4") +
       theme_classic() +
       scale_y_continuous(limits = c(0, 2500), expand = c(0.03, 0), breaks = seq(0, 2500, by=500)) +
       scale_x_continuous(expand = c(0.03, 0), breaks = c(seq(0, 3100, by = 500),722)) +
       theme(axis.title = element_text(family = "serif", color = "black", size = 12),
             axis.text  = element_text(family = "serif", color = "black", size = 12),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "white"),
             panel.grid.minor = element_line(colour = "white")) +
       xlab("eccDNA size (bp)") +ylab("Number of eccDNA")
p_2
ggsave(plot = p_2, filename = "Cancer_hsa_under3000bp.pdf", units = "cm", width = 18, height = 12.2)
dev.off()
#开始作图：智人的Cancer组别的特定检测方法的直方图
jpeg(filename = "Cancer_PredictedManually_hsa.jpg", res = 800, units = "cm", width = 9, height = 6)
p_3 <- ggplot(ha_m2) +
       geom_histogram(aes(x = Length), binwidth = 1, fill = "royalblue4") +
       theme_update() +
       scale_y_continuous(limits = c(0, 60), expand = c(0.01, 1), breaks = seq(0, 40, by = 10)) +
       scale_x_continuous(expand = c(0.03, 0.01), breaks = seq(0, 3000, by = 400)) +
       theme(axis.title = element_text(family = "serif", color = "black", size = 18),
             axis.text  = element_text(size = 16, family = "serif", color = "black"),
             plot.margin=unit(c(1,2,1,1),'lines')) +
       xlab("eccDNA size (bp)") +ylab("eccDNA counts")
p_3
dev.off()
#开始作图：小鼠的Healthy组别的直方图
jpeg(filename = "Healthy_mmu.jpg", res = 800, units = "cm", width = 9, height = 6)
p_4 <- ggplot(ma_dA) +
       geom_histogram(aes(x = Length), binwidth = 1, fill = "royalblue4") +
       #theme_update() +
       theme_classic() +
       scale_y_continuous(expand = c(0.02,0), breaks = seq(0, 10000,by = 1000)) +
       scale_x_continuous(expand = c(0.02,0), breaks = seq(0, 3000, by = 400)) +
       theme(axis.title = element_text(family = "serif", color = "black", size = 18),
             axis.text  = element_text(family = "serif", color = "black", size = 18)) +
       xlab("eccDNA size (bp)") +ylab("eccDNA counts")
p_4
dev.off()


#开始统计数据：查看Length分布的峰值
uniqT <- as.data.table(ha_m2)
uniqT <- uniqT[, .N, by=c("Length")]
setkey(uniqT,N)
uniqT[uniqT$N==max(uniqT$N),]





