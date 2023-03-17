rm(list=ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_region")
library("ggplot2")
library("ggsci")
ha <- read.csv("/disk5/shy/1shyJOB/7eccdna/raw/hsa2.txt",
               sep="\t",stringsAsFactors = F,header = T)
#View(colnames(ha))
ab <- cbind(ha[,grepl("frac",colnames(ha))],ha[,c(5,6,48)])
a_H <- ab[grepl("Healthy",ab$Disease),]
a_C <- ab[grepl("[^Healthy]",ab$Disease),]
#colnames(a_C)
a_region <- data.frame()
for (i in c("X5.UTR_frac","X3.UTR_frac",
            "exon_frac","intron_frac","Promoter_frac")) {
  a_region["Healthy", i] <- sum(a_H[,i])/length(rownames(a_H))
  a_region["Cancer", i] <- sum(a_C[,i])/length(rownames(a_C))
}
for (i in seq(2)) {
  print(i)
  a_region[i,"other"] <- 1-sum(a_region[i,1:5])
}
a_regionT <- as.data.frame(t(a_region))
rownames(a_regionT) <- c("5' UTR","3' UTR","Exon","Intron","Promoter","Intergenic region")

a_regionT[,c("region","HH","CC","q","w")] <- data.frame(
                           rownames(a_regionT),rep("Healthy",6),
                           rep("Cancer",6),seq(12)[1:6],seq(12)[7:12])
r_1 <- a_regionT[,c(4,3,1)]
rownames(r_1) <- a_regionT[,6]; colnames(r_1) <- seq(3)
r_2 <- a_regionT[,c(5,3,2)]
rownames(r_2) <- a_regionT[,7]; colnames(r_2) <- seq(3)
a_Plot <- rbind(r_1,r_2)
colnames(a_Plot) <- c("clas","regio","valu")
#作图
#pdf("/disk5/shy/1shyJOB/7eccdna/ecc_region/plotAnnoBar.pdf",width = 12,height = 5)
#jpeg("/disk5/shy/1shyJOB/7eccdna/ecc_region/plotAnnoBar.jpg",width = 24,height = 8,res = 600,units = "cm")
a_Plot$regio <- factor(a_Plot$regio, 
                       levels = c("Intergenic region","Intron","Exon","3' UTR","5' UTR","Promoter"))
p_region <- ggplot(a_Plot,aes(x=clas,y=valu,fill=regio)) +
            theme_classic() +
            geom_bar(width = 0.8,stat = "identity") +
            ylab("Fraction of Genomic region in eccDNA") +
            #scale_fill_discrete() +
            scale_fill_nejm() +
            theme(axis.text = element_text(family = "serif",
                                           size = 12,colour = "black"),
                  legend.text = element_text(family = "serif",
                                           size = 12,colour = "black"),
                  axis.title.x = element_text(family = "serif",
                                           size = 14,colour = "black"),
                  axis.title.y = element_blank(),
                  axis.line = element_line(size = 0.5),
                  legend.title = element_text(family = "serif",
                                           size = 12,colour = "black"),
                  legend.background = element_rect(fill = "white"),
                  panel.background = element_rect(fill = "white"),
                  panel.grid.major = element_line(colour = "white"),
                  panel.grid.minor = element_line(colour = "white")) +
            guides(fill=guide_legend(reverse = T,title = "Genomic region")) +
            coord_flip()
p_region
#ggsave("region.jpg",dpi = 800,units = "cm",width = 22,height = 6)
ggsave("region.pdf",units = "cm",width = 22,height = 6,plot = p_region)

dev.off()


as.numeric(a_H[2,][c("X5.UTR_frac","X3.UTR_frac","exon_frac","intron_frac","Promoter_frac")])
a_C$intergenic <- apply(a_C, 1, function(x){
  1 - sum(as.numeric(x[c("X5.UTR_frac","X3.UTR_frac",
              "exon_frac","intron_frac","Promoter_frac")]))
})




table_s1 <- data.frame()
for (i in c("X5.UTR_frac","X3.UTR_frac",
            "exon_frac","intron_frac","Promoter_frac","intergenic")) {
  table_s1["Healthy",i] <- round(sum(a_H[,i])/length(rownames(a_H)),5)
  table_s1["h_se",i]    <- round(sd(a_H[,i])/sqrt(length(rownames(a_H))),5)
  table_s1["Cancer",i]  <- round(sum(a_C[,i])/length(rownames(a_C)),5)
  table_s1["c_se",i] <- round(sd(a_C[,i])/sqrt(length(rownames(a_C))),5)
}
sum(table_s1["Healthy",])
