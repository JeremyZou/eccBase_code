#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/supply/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","tidyverse")
for (i in pks){library(i, character.only = T)}
rm(i,pks)

#开始导入数据:智人的repClass数据
rS_h <- read.csv("repClass_hsa.txt", header = F, stringsAsFactors = F, sep = "\t")
#开始导入数据:小鼠的repClass数据
rS_m <- read.csv("repClass_mmu.txt", header = F, stringsAsFactors = F, sep = "\t")


#开始作图：智人的repClass分布（条形图）
rS_h$V2 <- factor(rS_h$V2, levels = c(rev(as.character(rS_h$V2))[2],
                                      rev(as.character(rS_h$V2))[1],
                                      rev(as.character(rS_h$V2))[3:18]))
options(scipen=200)
p_1 <- ggplot(data = rS_h) +
      geom_bar(mapping = aes(x=V2,y=V1), fill="royalblue4",stat = "identity",width = 0.5) +
      theme_classic() +
      labs(x="", y="Number of eccDNAs",title = "") +
      scale_y_continuous(limits = c(0,400000), 
                         breaks = c(seq(0,400000,by=100000))) +
      scale_x_discrete(labels = c("Without_repeat", rev(as.character(rS_h$V2))[1], 
                                    rev(as.character(rS_h$V2))[3:18])) +
      theme(axis.title  = element_text(family = "serif", size = 12),
            axis.text   = element_text(family = "serif", size = 12, 
                                       angle = 30, hjust = 0.9, vjust = 0.8),
            plot.title  = element_text(family = "serif", size = 12, hjust = 0.5)) +
      geom_text(aes(x = V2, label = V1, y = V1 + 40000, family = "serif"), 
                angle = 90, size=4) 
p_1
getwd()
ggsave(p_1, filename = "repClass_hsa.jpg", device = "jpeg", dpi = 800, units = "cm", width = 28, height = 12)
ggsave(p_1, filename = "repClass_hsa_bar.pdf", device = "pdf", units = "cm", width = 28, height = 10)

#开始作图：小鼠的repClass分布（条形图）
rS_m$V2 <- factor(rS_m$V2, levels = c(rev(as.character(rS_m$V2))))
p_2 <- ggplot(data = rS_m) +
       geom_bar(mapping = aes(x=V2,y=V1), fill="royalblue4",stat = "identity", width = 0.5) +
       theme_classic() +
       labs(x="", y="Number of eccDNAs",title = "") +
       scale_y_continuous(limits = c(0,350000), 
                          breaks = c(seq(0,350000,by=100000))) +
       scale_x_discrete(labels = c("Without_repeat", 
                                   rev(as.character(rS_m$V2))[2:17])) +
       theme(axis.title  = element_text(family = "serif", size = 12),
             axis.text   = element_text(family = "serif", size = 12, 
                                        angle = 30, hjust = 0.9, vjust = 0.8),
             plot.title  = element_text(family = "serif", size = 12, hjust = 0.5)) +
       geom_text(aes(x = V2, label = V1, y = V1 + 30000, family = "serif"), 
                 angle = 90, size=4)
p_2
ggsave(p_2, filename = "repClass_mmu.jpg", device = "jpeg", dpi = 800, 
       units = "cm", width = 28, height = 12)
ggsave(p_2, filename = "repClass_mmu_bar.pdf", device = "pdf", units = "cm", width = 28, height = 10)


pieX <- data.frame(all_gF = c(0.4162,0.5838), protein_gF = c(0.4559,0.5441), 
                   nc_gF  = c(0.1458,0.8542), lnc_gF = c(0.1504,0.8496), 
                   pseudo_gF = c(0.0488,0.9512), nnam = c("hii","hii"),
                   filcol =c("gene","notG"))
gf1 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=all_gF,fill=filcol), stat = "identity") +
  scale_fill_manual(values = c("royalblue4", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf1
getwd()
ggsave(gf1, filename = "repClass_hsa_pie.pdf", device = "pdf", units = "cm", width = 4, height = 4)
ggsave(gf1, filename = "repClass_mmu_pie.pdf", device = "pdf", units = "cm", width = 4, height = 4)




