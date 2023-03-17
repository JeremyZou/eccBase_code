#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_repeat/")
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","sampling")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#
#
#开始整理数据：读取所有text文件，并进行统计
aa = "LTR"
aa = "LINE"
aa = "SINE"
okk <- assign(paste0("d_",aa),
       data.frame(dadsa <- seq(14), adaw <- seq(14), rewr <- seq(14)))
colnames(okk) <- c("nm", "all", "ltr")
f_all <- list.files(path = "./tempTissues_mmu",pattern = ".txt", full.names = T)
for (i in seq(length(list.files(path = "./tempTissues_mmu",pattern = ".txt", full.names = T)))) {
  temp       <- read.csv(f_all[i], header = F, sep = "\t", stringsAsFactors = F)
  okk$nm[i]  <- str_split(basename(f_all[i]), pattern = "\\.")[[1]][1]
  okk$all[i] <- length(temp[,1])
  okk$ltr[i] <- length(temp[grep(pattern = aa, x = temp[,25], perl = T),1])
}
okk$ltr_frac <- okk$ltr/okk$all*100
okk <- arrange(okk, desc(ltr_frac))
okk$group <- aa
data <- okk[,c(1,4,5)]

#开始作图：堆叠柱状图
data_LTR <- data
set.seed(111); xn <- sample(1:14,14)
data$nm <- factor(data$nm, levels = data_LTR$nm[xn])
p_ltr <- ggplot(data, aes(x=group, weight=ltr_frac, fill=nm)) +
         geom_bar(position = "dodge", show.legend = F, width = 1) +
         theme_classic() +
         scale_fill_d3(palette = "category20") +
         ylab(label = "eccDNAs (%)") +
         xlab(label = "repClass (LTR)") +
         #ggtitle(label = "Percentage of eccDNAs possess LTR across tissues of mmu") +
         theme(axis.text = element_text(family = "serif", size = 10, color = "black"),
               axis.text.x = element_blank(),
               axis.title = element_text(family = "serif", size = 10, colour = "black"),
               plot.title = element_text(family = "serif", size = 12, color = "black", hjust = 0.5),
               plot.margin = unit(c(0.2,0.2,0.2,0.4),"cm")) +
         scale_y_continuous(expand = c(0.04,0),limits = c(0,22),breaks = seq(0,22,5)) +
         scale_x_discrete(expand = c(0.53,0))
         #+geom_text(data2,family="serif", aes(x=group,y=frac+0.015,color=nm, label=c(round(frac,3)[1:3],"0.130",round(frac,3)[5:13],"0.100")),angle=90)
p_ltr

data$nm <- factor(data$nm, levels = data_LTR$nm[xn])
p_line <- ggplot(data, aes(x=group, weight=ltr_frac, fill=nm)) +
  geom_bar(position = "dodge", show.legend = F, width = 1) +
  theme_classic() +
  scale_fill_d3(palette = "category20") +
  ylab(label = "eccDNAs (%)") +
  xlab(label = "repClass (LINE)") +
  theme(axis.text = element_text(family = "serif", size = 10, color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_text(family = "serif", size = 10, colour = "black"),
        plot.title = element_text(family = "serif", size = 12, color = "black", hjust = 0.5)
        ) +
  scale_y_continuous(expand = c(0.04,0),limits = c(0,16),breaks = seq(0,16,5)) +
  scale_x_discrete(expand = c(0.53,0))
p_line

data$nm <- factor(data$nm, levels = data_LTR$nm[xn])
p_sine <- ggplot(data, aes(x=group, weight=ltr_frac, fill=nm)) +
          geom_bar(position = "dodge", show.legend = T, width = 1) +
          theme_classic() +
          scale_fill_d3(palette = "category20") +
          ylab(label = "eccDNAs (%)") +
          xlab(label = "repClass (SINE)") +
          labs(fill="Tissues in mmu") +
          theme(axis.text   = element_text(family = "serif", size = 10, color = "black"),
                axis.text.x = element_blank(),
                legend.key.size = unit(0.2,"cm"),
                legend.title = element_text(family = "serif", size = 10, color="black"),
                legend.text  = element_text(family = "serif", size = 10, color="black"),
                axis.title   = element_text(family = "serif", size = 10, colour = "black"),
                plot.title   = element_text(family = "serif", size = 12, color = "black", hjust = 0.5)
                ) +
          scale_y_continuous(expand = c(0.04,0),limits = c(0,35),breaks = seq(0,35,5)) +
          scale_x_discrete(expand = c(0.53,0))
p_sine

myplot <- p_ltr + p_line + p_sine + plot_layout(ncol = 3, widths = c(4,4,4))
myplot
#ggsave(myplot,filename = "./tempTissues_mmu/1_howManyecc.jpg", dpi = 800, units = "cm", width = 20, height = 6)
ggsave(myplot,filename = "./tempTissues_mmu/1_howManyecc.pdf", units = "cm", width = 20, height = 6)
getwd()

print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")
print("Over~Over~Over")





















#开始绘图：徒手绘制图例
ggplot() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white",colour = "white"),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  annotate("rect", xmin = 1, xmax = 2, ymin=1, ymax=2, fill="white") +
  annotate("rect", xmin = 1, xmax = 1.12, ymin=1.6, ymax=1.95, fill="indianred3") +
  annotate("rect", xmin = 1, xmax = 1.12, ymin=1.0, ymax=1.35, fill="gray60") +
  annotate("text", x = 1.52, y = 1.77, label="overlapping with LTR", family="serif",size=4.5) +
  annotate("text", x = 1.60, y = 1.24, label="non-overlapping with LTR", family="serif", size=4.5)
ggsave(filename = "./tempTissues_mmu/LTR_tissues_legend.jpg", dpi = 800, units = "cm", width = 6, height = 2)

#开始作图：按照示例的方法作图（堆叠柱状图）
data$group <- factor(data$group, levels = rev(c("ltr","noltr")))
ggplot(data = data, aes(x=nm, weight=frac, fill = group)) +
  geom_bar(position = "stack", show.legend = T, width = 0.9) +
  scale_fill_nejm()





