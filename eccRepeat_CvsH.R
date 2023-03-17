#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
scales::show_col(pal_nejm()(8))


#开始整理数据
ha <- read.csv("../raw/hsa2.txt",header = T,sep = "\t",stringsAsFactors = F)
hb_f <- cbind(ha[,grepl("frac",colnames(ha))],ha[,c(1,5,6,48)])
data <- hb_f[,c("repeat_frac","Disease")]
he   <- data %>% dplyr::filter(Disease == "Healthy"); he$Disease <- "Healthy"
ce   <- data %>% dplyr::filter(Disease != "Healthy"); ce$Disease <- "Cancer"
fp   <- rbind(he,ce)
#开始作图
#1、ggplot风格提琴图
fp$Disease <- factor(fp$Disease, levels = c("Healthy", "Cancer"))
jpeg(filename = "repeat_hsa_HvsC.jpg", res = 800,units = "cm", width = 12, height = 10)
p_1 <- ggplot(data = fp) + 
       theme_classic2() +
       geom_violin(aes(x = Disease, y = repeat_frac),
                   width = 0.9, fill = "red",color="white") +
       stat_compare_means(aes(x=Disease,y=repeat_frac)) +
       scale_y_continuous(limits = c(0, 1))
p_1


ggviolin(data = fp, x = "Disease", y = "repeat_frac")







#2、vioplot风格提琴图
jpeg(filename = "repeat_hsa_HvsC.jpg", res = 800,units = "cm", width = 9, height = 7)
pdf(file = "./CvsH/repeat_hsa_HvsC.pdf",width = 9, height = 7,family = "serif")
par(family="serif", ps = 12, col=c("white"), mar=c(2.2,2.2,0,0))
vioplot::vioplot(fp[which(fp$Disease=="Healthy"),1],fp[which(fp$Disease=="Cancer"),1],
                 names = c("Healthy","Cancer"),
                 axes  = F, 
                 col = c("indianred3"), cex.axis = 1,
                 family = "serif", ylim = c(0,1.15), side=c(1,2))
dev.off()


#
#
#
#
#
#
#
#
#统计
#统计学计算：平均值、标准差、中位数、维尔康克森秩和检验显著性差异
print("平均数");mean(he[,1]);print("标准差");sd(he[,1]);print("中位数");median(he[,1])
print("平均数");mean(ce[,1]);print("标准差");sd(ce[,1]);print("中位数");median(ce[,1])
#Healthy组别中重复序列占比超过80%的eccDNA达到了43.34%
length(he[he[,1]>0.8,1])/145399;length(he[he[,1]>0.5,1])/145399
#Cancer 组别中重复序列占比超过80%的eccDNA达到了26.27%
length(ce[ce[,1]>0.8,1])/591436;length(ce[ce[,1]>0.5,1])/591436
compare_means(data = fp, formula = repeat_frac~Disease, method = "wilcox.test",
              paired = F,p.adjust.method = "fdr")
wilcox.test(fp[which(fp$Disease=="Healthy"),1],fp[which(fp$Disease=="Cancer"),1],
            alternative = "greater")
