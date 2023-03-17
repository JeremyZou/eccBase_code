#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_gc/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap")
for (i in pks){library(i, character.only = T)}
rm(i,pks)


#开始整理数据：智人的Healthy、Cancer组别的GC含量
ha <- read.csv("../raw/hsa2.txt",header = T,sep = "\t",stringsAsFactors = F)
gc <- ha[,c("GC_frac","Disease")]
aH         <- gc[grepl("Healthy",gc$Disease),]
aH$Disease <- "Healthy"
aC         <- gc[grep(pattern = "[^Healthy]",x = gc$Disease),]
aC$Disease <- "Cancer"
#开始整理数据：合并Healthy、Cancer组别的GC含量数据
pData <- rbind(aH,aC)
#开始整理数据：来自ClinVar、ESP6500、1000 Genome Project的GC含量数据
cl <- read.csv("clinvar_delBenign_pop_5bp.gc",
               header = T, stringsAsFactors = F, sep = "\t")
es <- read.csv("ESP6500_del_5bp.gc",
               header = T, stringsAsFactors = F, sep = "\t")
gp <- read.csv("1000GP_5bp.gc",
               header = T, stringsAsFactors = F, sep = "\t")
cl <- cl[,c(5,1)]; cl$X.1_usercol <- "ClinVar"; colnames(cl) <- c("GC_frac","Disease")
es <- es[,c(5,1)]; es$X.1_usercol <- "ESP6500"; colnames(es) <- c("GC_frac","Disease")
gp <- gp[,c(5,1)]; gp$X.1_usercol <- "1000 GP"; colnames(gp) <- c("GC_frac","Disease")
#开始整理数据：合并Healthy、Cancer、ClinVar、ESP6500、1000 Genome Project组别的GC含量数据
pData <- rbind(aH,aC,cl,es,gp)


#开始作图：ggplot风格提琴图 + 箱线图【所有组别】
pData$Disease <- factor(pData$Disease, levels = c("Healthy","Cancer", "ClinVar","ESP6500","1000 GP"))
pData$Disease <- factor(pData$Disease, levels = c("Healthy","Cancer"))
p_1 <- ggplot(data = pData, aes(x = Disease, y = GC_frac)) +
       geom_violin(fill="red") #+
       #geom_boxplot(width=0.1, outlier.shape = NA, fill="royalblue4", color="white")
p_1
#细节修改：使用geom_segment以覆盖的形式更改geom_boxplot()的middle line的颜色
dat <- ggplot_build(p_1)$data[[1]]
p_1 +  geom_segment(data = dat, aes(x = xmin, xend = xmax, y = middle, yend = middle), 
                    color = "red", inherit.aes = F,size=1.5)

#开始作图：vioplot风格提琴图
jpeg(filename = "GC_hsa_HvsC.jpg", res = 800,units = "cm", width = 9, height = 7)
pdf(file = "GC_hsa_HvC.pdf",width = 9, height = 7)
par(family="serif", ps = 12, col=c("white"), mar=c(2.2,2.2,0,0))
vioplot::vioplot(aH$GC_frac,aC$GC_frac,
         names = c("Healthy","Cancer"),
         axes  = F, 
         col = c("indianred3"), cex.axis = 1,
         family = "serif", ylim = c(0,1.15), side=c(1,2),
         ylab="sda")
dev.off()


























ggplot(aC,aes(GC_frac))+geom_density() +
  geom_density(data=aH,aes(GC_frac),color="red") +
  geom_density(data=es,aes(GC_frac),color="blue") +
  geom_density(data=gp,aes(GC_frac),color="cyan")
  
p <- vioplot(data=pData,GC_frac~Disease,
             col = c("gold"),cex.names = 1,
             family = "serif",ylim = c(0,1.1),line = c(0.1,0.9,0.5,0.5))

vioplot(aH$GC_frac)
median(aH$GC_frac)
summary(aC$GC_frac)
summary(cl$X8_pct_gc)
summary(es$X8_pct_gc)
summary(gp$X8_pct_gc)
#wilcox.test(aH$GC_frac,aC$GC_frac,)



data_one <- rnorm(100)
data_two <- rnorm(50, 1, 2)
vioplot(data_one, data_two, col=c("red","blue"), names=c("data one", "data two"),
        main="data violin", xlab="data class", ylab="data read")
legend("topleft", fill=c("red","blue"), legend=c("data one", "data two"))






