rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_disease/")
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggsunburst)
library(data.table)
library(ggnewscale)
library(gridExtra)
##################
rawF <- read.csv("../raw/mmu.txt",sep = "\t",stringsAsFactors = F,header = T)
aa <- rawF[,c("Disease","Cell.Tissues","Identify_method")]
ad <- read.csv("dtcm_mmu",sep = "\t",header = F,stringsAsFactors = F)
df <- data.frame()
ti = 0
for (line in rownames(ad)){
  ti = ti+1
  print (ti)
  cc=paste0("row_",ti)
  print(cc)
  cc <- aa[grep(pattern = ad[line,]$V1, aa$Disease),]
  cc$Disease <- ad[line,]$V1
  cc <- cc[grep(pattern = ad[line,]$V2, cc$Cell.Tissues),]
  cc$Cell.Tissues <- ad[line,]$V2
  cc <- cc[which(cc$Identify_method == ad[line,]$V3),]
  cc$Identify_method <- ad[line,]$V3
  df <- rbind(df, cc)
}
dfd <- as.data.frame(table(df$Disease))
dfc <- as.data.frame(table(df$Cell.Tissues))
dfm <- as.data.frame(table(df$Identify_method))
df01 <- df
ti=0
for(kin in rownames(dfc)){
  ti=ti+1; print(ti)
  df01$Cell.Tissues <- gsub(dfc$Var1[ti],dfc$Freq[ti],df01$Cell.Tissues)
}
df_f<-cbind(df,df01)
colnames(df_f) <- c("disease","celltissue","imeth","couD","couC","couM")
write.table(df_f,"df_f",col.names = T,row.names = F,quote = F,sep = "\t")
df_cell <- df_f

#准备作图的数据
Df_f <- as.data.table(df_cell)
sb_0 <- Df_f[, .(xmin = 0.0, xmax = 1.2, ymin = 0, ymax = .N)]
sb_1 <- Df_f[, .N, by = c("celltissue")]
#排序方法一
#sb_1$celltissue <- factor(sb_1$celltissue, levels = objord)
#sb_1 <- sb_1[order(sb_1$celltissue),]
#排序方法二
setkey(sb_1, N)
sb_1[, `:=`(xmin = 1.2, xmax = 2.5, ymin = cumsum(shift(N, fill = 0)), ymax = cumsum(N))]
#digit <- sprintf("%0.4f",sb_1$N/sum(sb_1$N))
sb_1$frac <- paste(round(sb_1$N/sum(sb_1$N),4)*100,"%")
#################
#start to plot
g_rec <-
  ggplot() +
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax) +
  geom_rect(data = sb_0, fill = "white") +
  geom_rect(data = sb_1, mapping = aes(fill = celltissue), colour="white") +
  scale_fill_manual(breaks = sb_1$celltissue,
                    values = rev(c(pal_d3(palette = "category20")(20)[1:14])),
                    labels = c(paste(sb_1$celltissue," (",sb_1$frac,")"))) +
  #  geom_text_repel(data = sb_1, 
  #                  aes(x=3.5,y=(ymax-ymin)/2,label=frac,hjust=0,
  #                                   colour="black",family="serif"),size=4) +
  guides(fill = guide_legend(order = 1, nrow = 14, title = "Tissue & cell line",reverse = T))
g_rec
g_sb_04 <-
  g_rec +
  coord_polar(theta = "y",direction = -1) +
  ggtitle(label = "         mmu eccDNA distribution by tissue & cell line") +
  theme_classic() +
  theme(legend.position = c(1.28,0.46),
        legend.text = element_text(family = "serif", size = 16),
        legend.title = element_text(family = "serif", size = 16),
        plot.margin = unit(c(0,10,0,0),"cm"),
        plot.title = element_text(family = "serif", size = 18,
                                  hjust = 0.5)) +
  theme(axis.title = element_blank(),
        axis.text  = element_blank(),
        axis.ticks = element_blank(),
        axis.line  = element_blank()) +
  annotate("text", label = "mmu (481831)\n---------------\nIn total / Healthy", 
           x=0, y=0.5, family="serif", size=6)
g_sb_04
#ggsave("rectPie_cell_mmu.jpg",dpi=800)
ggsave(g_sb_04,filename = "AI/rectPie_cell_mmu.pdf",units = "cm", width = 24, height = 15)



#组合成图2、、、、、失败了
myplot_ltr <- g_sb_01 + g_sb_02 + plot_layout(ncol = 2, widths = c(20,20))
myplot_ltr
getwd()
ggsave(myplot_ltr,filename = "AI/test.pdf",units = "cm", width = 30, height = 30)













