#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_gene/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","ggupset","tidyverse",
         "devtools","easyGgplot2","ggforce","UpSetR")
for (i in pks){library(i, character.only = T)}
rm(i,pks)

#开始整理数据：智人的基因Count的相关数据
ha <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/hsa2.txt", sep="\t", stringsAsFactors = F, header = T)
hb <- cbind(ha[,grepl("count",colnames(ha))],ha[,c(5,6,48)])
#临时
hb  <- cbind(ha[,grepl("frac",colnames(ha))],ha[,c(5,6,48)])
h_H <- hb %>% dplyr::filter(Disease == "Healthy")
h_C <- hb %>% dplyr::filter(Disease != "Healthy")
median(h_H$repeat_frac)
median(h_C$repeat_frac)

uniqA <- as.data.frame(apply((as.data.frame(table(hb$allGene_count))),2,function(x){
  as.integer(x)
}))
uniqA$Lo <-log2(uniqA$Freq +1)
#开始整理数据：智人的基因Count的相关数据（按照疾病分组）
h_H <- hb %>% dplyr::filter(Disease == "Healthy")
h_C <- hb %>% dplyr::filter(Disease != "Healthy")
uniqH <- as.data.frame(apply((as.data.frame(table(h_H$allGene_count))),2,function(x){
  as.integer(x)
}))
uniqH$Lo <-log2(uniqH$Freq +1)
uniqC <- as.data.frame(apply((as.data.frame(table(h_C$allGene_count))),2,function(x){
  as.integer(x)
}))
uniqC$Lo <-log2(uniqC$Freq +1)

#开始作图：智人的基因Count的柱状图
p_1 <- ggplot(data = uniqA) +
       theme_classic() +
       ggtitle(label = "Frequency distribution of hsa \neccDNA and located gene") +
       theme(axis.title = element_text(family = "serif",size = 20),
             axis.text  = element_text(family = "serif",size = 20),
             plot.title = element_text(family = "serif",size = 20,hjust = 0.5),
             axis.line  = element_line(colour = "black",size = 0.8),
             axis.ticks = element_line(colour = "black",size = 0.8)) +
       geom_bar(aes(x=Var1,y=Lo), stat = "identity",width = 0.2,color="indianred3") +
       labs(y="Log2(Number of\n eccNDAs + 1)",x="Number of eccDNA genes") +
       scale_y_continuous(limits = c(-0.2,20),breaks = seq(0,20,by=4)) +
       scale_x_continuous(limits = c(-0.2,4200),breaks = c(seq(0,3500,by=500),4041),
                          expand = c(0.015,0))
p_1
ggsave(p_1, filename = "allGene_hsa.jpg",units = "cm",width = 16,height = 10, dpi = 800)
#开始作图：智人的基因Count的直方图 + 频率折线图
jpeg(filename = "allGene_histogram_bar_hsa.jpg",res = 800,units = "cm", width = 8, height = 4)
p_2 <- ggplot() +
       theme_classic() +
       geom_bar(data = uniqA, aes(x=Var1,y=Lo), stat = "identity",width = 0.2,color="indianred3",alpha=1) +
       #coord_cartesian(xlim = c(0,2000)) +
       scale_y_continuous(expand = c(0.02,0),limits = c(-0.2,20),breaks = seq(0,20,by=4),
                          sec.axis = sec_axis(trans = ~.,breaks = seq(0,20,by=4),
                                              labels = c(0,paste0(seq(0,20,by=4)[2:6]*3,"%")),
                                              name = "Frequency")) +
       scale_x_continuous(limits = c(),breaks = c(seq(0,3500,by=1000),4041),
                         expand = c(0.02,0)) +
       geom_freqpoly(data = hb[,c(1,2)], aes(allGene_count, y = ..count../736775*100/3), 
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
       xlab("Number of eccDNA genes") +ylab("Log2(Number of eccDNAs + 1)")
p_2
ggsave(plot = p_2, filename = "allGene_histogram_bar_hsa_v2.pdf", units="cm", width = 8, height = 6)
dev.off()


#uniqA <- as.data.table(uniqA)
#setkey(uniqA,Var1)
#1-uniqA[which(uniqA$Var1 == 0),2]/sum(uniqA[,2])
pieX <- data.frame(all_gF = c(0.6354,0.3646), protein_gF = c(0.4559,0.5441), 
                   nc_gF  = c(0.1458,0.8542), lnc_gF = c(0.1504,0.8496), 
                   pseudo_gF = c(0.0488,0.9512), nnam = c("hii","hii"),
                   filcol =c("gene","notG"))
gf1 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=all_gF,fill=filcol), stat = "identity") +
                 scale_fill_manual(values = c("royalblue4", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf1
getwd()
ggsave(plot = gf1, filename = "./supply/pie_hsa_gene.pdf", units = "cm", width = 4, height = 4)

gf2 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=protein_gF,fill=filcol), stat = "identity") +
                 scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf3 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=nc_gF,fill=filcol), stat = "identity") +
                 scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf4 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=lnc_gF,fill=filcol), stat = "identity") +
                 scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
gf5 <- ggplot() +geom_bar(data = pieX, mapping = aes(x=nnam,y=pseudo_gF,fill=filcol), stat = "identity") +
                 scale_fill_manual(values = c("indianred", "gray90")) +theme(legend.position = "none") +coord_polar(theta = "y",direction = 1) +theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),axis.title = element_blank(),panel.background = element_rect(fill = "transparent"))
#gf6 <- ggplot() +geom_rect()
jpeg(filename = "diffGene_frac_hsa_1.jpg", width =32, units = "cm", height = 4, res = 800)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 5)))
print(gf1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(gf2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(gf3, vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(gf4, vp = viewport(layout.pos.row = 1, layout.pos.col = 4))
print(gf5, vp = viewport(layout.pos.row = 1, layout.pos.col = 5))
#ggsave(filename = "diffGene_frac_hsa.jpg", dpi=800)
dev.off()
#ggsave(gf1,filename = "allGene_06354_hsa.jpg", dpi=800)
#ggsave(gf,filename = "protein_04559_hsa.jpg", dpi=800)
#ggsave(gf,filename = "nc_01458_hsa.jpg", dpi=800)
#ggsave(gf,filename = "lnc_01504_hsa.jpg",dpi=800)
#ggsave(gf,filename = "pseudo_00488_hsa.jpg",dpi=800)



#temp <- ggplot() +geom_bar(aes(x=c("sa"),y=c(0.7,0.3),fill=c("sa","er")),stat = "identity") +
#  coord_polar(theta = "y",direction = -1.5) +
#  theme(axis.line = element_blank(),axis.ticks = element_blank(),
#        axis.text = element_blank(),axis.title = element_blank(),
#        panel.background = element_rect(fill = "transparent")) +
#  scale_fill_manual(values=c("white","indianred")) +
#  guides(fill="none")
#ggsave(temp,filename = "temp.jpg",device = "jpeg",dpi=600)




#uniqA[which(uniqA[,"Var1"]>4000),]
#uniqA[which(uniqA[,"Var1"]>1950 & uniqA[,"Var1"]<2050),]
#length(hb[which(hb$protein_count > 0),][,1])
#length(hb[which(hb$ncRNA_count > 0),][,1])
#length(hb[which(hb$lncRNA_count > 0),][,1])
#length(hb[which(hb$pseudos_count > 0),][,1])
#ha$ncRNA_name[seq(10)]




#a_H <- hb[grepl("Healthy",hb$Disease),]
#a_C <- hb[grepl("[^Healthy]",hb$Disease),]
#a_H$Disease <- "Healthy";a_C$Disease <- "Cancer"
#sd(a_H$X5.UTR_frac);mean(a_H$X5.UTR_frac)
#length(a_H[which(a_H$tss_count==0),][,"tss_count"])/length(a_H[,1])
#length(a_H[which(a_H$tss_count>0 & a_H$tss_count<10),][,"tss_count"])/length(a_H[,1])
#length(a_H[which(a_H$tss_count>=10),][,"tss_count"])/length(a_H[,1])
#length(a_C[which(a_C$tss_count==0),][,"tss_count"])/length(a_C[,1])
#length(a_C[which(a_C$tss_count>0 & a_C$tss_count<10),][,"tss_count"])/length(a_C[,1])
#length(a_C[which(a_C$tss_count>=10),][,"tss_count"])/length(a_C[,1])


#ggplot(a_C[,c("tss_count","Disease")],aes(x=tss_count)) +
#  geom_freqpoly(binwidth = 1)+
#  geom_histogram(binwidth = 1,fill="royalblue4") +
#  theme_update()

#ggplot(diamonds, aes(price, fill = cut)) +
#  geom_histogram(binwidth = 500)
#ggplot(diamonds, aes(price, colour = cut)) +
#  geom_freqpoly(binwidth = 500)


#da[,length(colnames(da))+1] <- apply(da[,c(1,2)], 1, function(x) {
#  #inds <- match(x[1], unicq$Var1)
#  inds <- match(as.character(x[1]), as.character(unicq$Var1))
#  ifelse(is.na(inds), x[1], unicq$Freq[inds]) 
#})


#.libPaths("/home/zyw/R/x86_64-pc-linux-gnu-library/3.5")
#download.file("https://github.com/kassambara/ggpubr/archive/refs/tags/v0.4.0.tar.gz",destfile = "./ggpur.tar.gz")
