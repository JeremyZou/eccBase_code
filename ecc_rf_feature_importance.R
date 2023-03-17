#开始组建环境
#rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_predict/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#开始整理数据：导入数据
rf_fea <- read.csv(file = "rf_feature_importance_XGB.txt", header = F, sep = "\t", stringsAsFactors = F, row.names = 1)
rf_fea$nm <- rownames(rf_fea)
rf_fea$we <- "qqqq"
rf_fea <- rf_fea[order(rf_fea$V2,decreasing = T),]
rf_fea$uu <- paste0(rf_fea$nm," (", round(rf_fea$V2,3),")")
rf_fea$uu[6:16] <- rf_fea$nm[6:16]
rf_fea$uu <- factor(rf_fea$uu, levels = rf_fea$uu)
p_2 <- ggplot(data = rf_fea) +
       geom_bar(mapping = aes(x=we,y=V2,fill=uu),stat = "identity", 
                position = "stack", width = 1) +
       labs(fill="Features") +
       ylab(label = "Importance score") +
       ggtitle(label = "The XGBoost feature importance") +
       theme_classic() +
       theme(plot.margin = unit(c(0.2,0.8,0.2,0.2),"cm"), 
             axis.text = element_text(family = "serif", size = 12, color = "black"),
             axis.title.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.text.x = element_text(colour = "white"),
             axis.text.y = element_text(family = "sans", size = 8, colour = "black"),
             axis.title.y = element_text(family = "sans", size = 8, color = "black"),
             plot.background = element_rect(fill="white",colour = "white"),
             plot.title = element_text(family = "sans", size = 8,  color = "black", hjust = 0),
             legend.text = element_text(family = "sans", size = 8, color = "black"),
             legend.title = element_text(family = "sans", size = 8,  color = "black")) +
       scale_y_continuous(limits = c(0,1.02),expand = c(0.02,0),breaks = seq(0,1,0.2)) +
       scale_fill_d3(palette = "category20")
p_2
ggsave(plot = p_1, filename = "rf_feature_importance.jpg", device = "jpeg", dpi = 800, units = "cm", width = 9.5, height = 11.5)

#合并AUC柱状图和特征重要性累积图
myplot_ltr <- p_1 + p_2 + plot_layout(ncol = 2, widths = c(8,2))
myplot_ltr
getwd()
ggsave(plot = myplot_ltr, filename = "predictUnionFig.jpg", device = "jpeg", 
       dpi = 800, units = "cm", width = 22, height = 12)

ggsave(plot = myplot_ltr, filename = "predictUnionFig.pdf", device = "pdf", 
       units = "cm", width = 22, height = 12)


