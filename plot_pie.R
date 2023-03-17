rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_disease/")
library(dplyr)
library(ggplot2)
#BiocManager::install("rlist")
library(rlist)
library(scales)
library(ggsci)
library(ggsunburst)
library(data.table)
library(ggnewscale)
library(gridExtra)
rawHF    <- read.csv("../raw/hsa2.txt",sep = "\t",stringsAsFactors = F,header = T)
rawMF    <- read.csv("../raw/mmu.txt",sep = "\t",stringsAsFactors = F,header = T)
hsa_dtcm <- read.csv("dtcm",sep = "\t",header = F,stringsAsFactors = F)
mmu_dtcm <- read.csv("dtcm_mmu",sep = "\t",header = F,stringsAsFactors = F)
diseN    <- unique(hsa_dtcm[,1])
cellN    <- unique(hsa_dtcm[,2])[grepl(unique(hsa_dtcm[,2]),pattern = "cell lines")]
tissueN  <- unique(hsa_dtcm[,2])[grep(unique(hsa_dtcm[,2]),pattern = "cell lines",invert=T)]
celmmuN  <- unique(mmu_dtcm[,2])
da       <- list()
ta       <- list()
for (i in diseN) {
  Le <- length(rawHF[grepl(rawHF$Disease,pattern = i),][,1])
  da <- rlist::list.append(da,i)
  ta <- rlist::list.append(ta,Le)
}
da   <- unlist(da)
ta   <- unlist(ta)
data <- data.frame(disName = da, disFreq = ta); data$ewai <- "zuo"
data[,1] <- as.character(data[,1])
data["88",] <- c("other",sum(data[which(data[,2]<3000),][,2]),"zuo")
data[,2] <- as.numeric(data[,2])
data <- data[which(data$disFreq>3000),]
data["99",] <- c("Inputt", sum(data[,2]), "koo")
data[,2] <- as.numeric(data[,2])
data
data <- data[order(data[,2], decreasing = T),]
data$disName <- factor(data$disName,levels = data$disName)
data$ewai <- factor(data$ewai, levels = c("koo","zuo"))


ggplot(data = data) +
  geom_bar(mapping = aes(x = ewai, y = disFreq, fill = disName),
           width=1.8,
           stat = "identity",position = "stack") +
  #scale_fill_d3(palette="category20") +
  #scale_fill_manual(values = c("gray93",colors()[c(63:81)])) +
  scale_fill_manual(values = c("red",pal_d3(palette = "category20")(19))) +
  coord_polar(theta = "y",direction = -1) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "gold",color = "black"),
        axis.title = element_blank(),axis.ticks = element_blank(),
        axis.text = element_blank(),axis.line = element_blank())



rainbow()


#Hi, here are the colors related
colors()
show_col(colors(distinct = T))
library("scales")
??show_col
show_col(pal_aaas(palette = "default")(5))
show_col(pal_igv()(51))
pal_aaas(palette = "default")(5)
pal_d3(palette = "category20")(20)
scales::show_col()
BiocManager::install("paletteer")
library("paletteer")

View(paletteer_packages)
