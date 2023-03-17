#开始组建环境
#rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_predict/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入数据
aucRE <- read.csv(file = "train_10fold_AUC_addXGB.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
aucRE$me <- apply(aucRE,1,function(x){mean(x)})
aucRE$sd <- apply(aucRE,1,function(x){sd(x)/sqrt(10)})

#开始作图
aucRE$wq <- rownames(aucRE)
aucRE$wq <- factor(aucRE$wq, levels = c("XGB","RF","MLP","SVC","LR","Bayes"))
p_1 <- ggplot(data = aucRE) +
       geom_bar(mapping = aes(x=wq,y=me),stat = "identity",fill=pal_jama()(6), width = 0.5) +
       geom_errorbar(aes(x=wq,y=me,ymin=me-sd,ymax=me+sd),width=0.2) +
       theme_classic() +
       ylab(label = "AUC") +
       ggtitle(label = "The mean AUC with 10-fold cross-validation") +
       theme(axis.title.x = element_blank(),
             axis.title = element_text(family = "sans", size = 8, color = "black"),
             axis.text = element_text(family = "sans", size = 8, color = "black"),
             plot.title = element_text(family = "sans", size = 8, color = "black", hjust = 0.5)) +
       scale_y_continuous(limits = c(0,1), breaks = seq(0,1,by = 0.2),expand = c(0.02,0)) +
       geom_text(size=3.2,aes(x = wq, label = c(0.946,"0.940",0.886,0.891,0.843,0.823), 
                     y = me + 0.05, family = "sans"))
p_1




ggsave(plot = p_1, filename = "auc_bar_XGB.jpg", device = "jpeg", dpi = 800, units = "cm", width = 11, height = 11.5)


