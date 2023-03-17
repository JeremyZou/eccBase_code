rm(list=ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_gene")
BiocManager::install("gg.gap")
library("ggplot2")
library("ggsci")
library("gg.gap")
library("scales")
library("easyGgplot2")
mm <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/mmu.txt",
               sep="\t",stringsAsFactors = F,header = T)
mb <- cbind(mm[,grepl("count",colnames(mm))],mm[,c(1,5,6,48)])
log2(length(mb[which(mb$repeat_count!=0),1]))
200547/481831




#1、for using geom_bar() to plot
uniqM <- as.data.frame(apply((as.data.frame(table(mb$allGene_count))),2,function(x){
  as.integer(x)
}))
uniqM$Lo <-log2(uniqM$Freq +1)
#2、for using geom_histogram() to plot
ab_pm <- mb[,c("allGene_count","eccDNA_ID")]
#1、开始作图:bar：柱状图
p_m <- ggplot(data = uniqM) +
       theme_classic() +
       ggtitle(label = "Frequency distribution of mmu eccDNA and located gene") +
       theme(axis.title = element_text(family = "serif",size = 16),
             axis.text = element_text(family = "serif",size = 16),
             plot.title = element_text(family = "serif",size = 16,hjust = 0.5)) +
       geom_bar(aes(x=Var1,y=Lo), stat = "identity",width = 0.2,color="indianred3",fill="indianred") +
       labs(y="Log2(Number of eccNDAs + 1)",x="Number of eccDNA genes") +
       scale_y_continuous(limits = c(-0.2,20),breaks = seq(0,20,by=2)) +
       scale_x_continuous(limits = c(-0.2,27),breaks = seq(0,27,by=2),expand = c(0.015,0))
p_m
ggsave(p_m, filename = "allGene_mmu.jpg",units = "cm",width = 8,height = 4, dpi = 800)
#
#2、开始作图:histogram：直方图 + 频率折线图
jpeg(filename = "allGene_histogram_bar_mmu.jpg",res = 800,units = "cm", width = 8, height = 4.1)
p_2m <- ggplot() +
  theme_classic() +
  geom_bar(data = uniqM, aes(x=Var1,y=Lo), stat = "identity",width = 0.2,
           color="indianred3",fill="indianred3",alpha=1) +
  #coord_cartesian(xlim = c(0,2000)) +
  scale_y_continuous(expand = c(0.02,0),limits = c(-0.2,20),breaks = seq(0,20,by=4),
                     sec.axis = sec_axis(trans = ~.,breaks = seq(0,20,by=4),
                                         labels = c(0,paste0(seq(0,20,by=4)[2:6]*3,"%")),
                                         name = "Frequency")) +
  scale_x_continuous(breaks = c(seq(0,28,by=4)),expand = c(-0.01,0)) +
  geom_freqpoly(data = ab_pm, aes(allGene_count, y = ..count../481831*100/3), 
                binwidth = 1, color="#1F77B4FF", size=1, alpha=1) +
  theme(axis.title.x       = element_text(family = "serif",color = "black",size = 10),
        axis.title.y.left  = element_text(family = "serif",color = "indianred3",size = 10),
        axis.title.y.right = element_text(family = "serif",color = "#1F77B4FF",size  = 10,
                                          angle = 90,vjust = 1),
        axis.text.y.left   = element_text(size = 10,family = "serif",color = "indianred3"),
        axis.text.y.right  = element_text(size = 10,family = "serif",color = "#1F77B4FF"),
        axis.line.y.left   = element_line(colour = "indianred3",size = 0.8),
        axis.line.y.right  = element_line(colour = "#1F77B4FF",size  = 0.8),
        axis.ticks.y.left  = element_line(colour = "indianred3",size = 0.8),
        axis.ticks.y.right = element_line(colour = "#1F77B4FF",size  = 0.8)) +
  xlab(expression(Number~of~eccDNA~genes)) +ylab("Log2(Number of eccDNAs + 1)")
p_2m
ggsave(p_2m, filename = "allGene_histogram_bar_mmu_v2.pdf",units = "cm", width = 8, height = 6.1)
dev.off()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
uniqM <- as.data.table(uniqM)
setkey(uniqM,Var1)
1-uniqM[which(uniqM$Var1 == 0),2]/sum(uniqM[,2]); uniqM[which(uniqM$Var1 == 0),2]/sum(uniqM[,2])
pieXm <- data.frame(all_gF = c(0.6035,0.3965), protein_gF = c(0.5496,0.4504), 
                   nc_gF  = c(0.0649,0.9351), lnc_gF = c(0.0686,0.9314), 
                   pseudo_gF = c(0.0057,0.9943), nnam = c("hii","hii"),
                   filcol =c("gene","notG"))
gf1 <- ggplot() +geom_bar(data = pieXm, mapping = aes(x=nnam,y=all_gF,fill=filcol), stat = "identity") +
  scale_fill_manual(values = c("royalblue4", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
ggsave(plot = gf1, filename = "./supply/pie_mmu_gene.pdf", units = "cm", width = 4, height = 4)


gf2 <- ggplot() +geom_bar(data = pieXm, mapping = aes(x=nnam,y=protein_gF,fill=filcol), stat = "identity") +
  scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf3 <- ggplot() +geom_bar(data = pieXm, mapping = aes(x=nnam,y=nc_gF,fill=filcol), stat = "identity") +
  scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf4 <- ggplot() +geom_bar(data = pieXm, mapping = aes(x=nnam,y=lnc_gF,fill=filcol), stat = "identity") +
  scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf5 <- ggplot() +geom_bar(data = pieXm, mapping = aes(x=nnam,y=pseudo_gF,fill=filcol), stat = "identity") +
  scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
#gf6 <- ggplot() +geom_rect()
jpeg(filename = "diffGene_frac_mmu_1.jpg", width =50, units = "cm", height = 12, res = 800)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 5)))
print(gf1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(gf2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gf3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(gf4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(gf5, vp = viewport(layout.pos.row = 1, layout.pos.col = 5))
#ggsave(filename = "diffGene_frac_hsa.jpg", dpi=800)
dev.off()

