#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/supply/");getwd()
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","ggupset","tidyverse",
         "devtools","easyGgplot2","ggforce","UpSetR")
for (i in pks){library(i, character.only = T)}
rm(i,pks)


1-length(hb[which(hb[,6]>0),6])/length(hb[,6])
1-length(mb[which(mb[,6]>0),6])/length(mb[,6])





#开始整理数据
pieX <- data.frame(all_gF     = c(0.6035,0.3965), #智人：0.6354,0.3646
                   all_reP    = c(0.4162,0.5838), #小鼠：0.7008,0.2992
                   protein_gF = c(0.4559,0.5441), 
                   nc_gF      = c(0.1458,0.8542), 
                   lnc_gF     = c(0.1504,0.8496), 
                   pseudo_gF  = c(0.0488,0.9512), 
                   nnam = c("hii","hii"), filcol =c("gene","notG"))
#开始作图
jpeg(filename = "allrep_frac_supply_mmu.jpg", width =8, units = "cm", height = 8, res = 800)
gf1 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=all_reP,fill=filcol), stat = "identity") +
       scale_fill_manual(values = c("black", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf1
dev.off()









#开始排版：一页多图
#grid.newpage()
#pushViewport(viewport(layout = grid.layout(1, 5)))
#print(gf1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
#print(gf2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
#print(gf3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
#print(gf4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
#print(gf5, vp = viewport(layout.pos.row = 1, layout.pos.col = 5))






