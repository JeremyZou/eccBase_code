#开始组建环境
rm(list = ls())
setwd("/disk5/shy/1shyJOB/7eccdna/ecc_gene/enrich/"); getwd(); list.files()
pks <- c("ggplot2","ggpubr","ggsci","scales","gg.gap","clusterProfiler",
         "topGO","Rgraphviz","pathview","org.Hs.eg.db","sampling")
for (i in pks){library(i, character.only = T)}
rm(i,pks)
#ls("package:org.Hs.eg.db")

#开始整理数据
i = 500
i = 1000
#开始整理数据：Cancer组别的基因
h_Cg   <- read.csv(file = "cancer_sortGene.rnk", stringsAsFactors = F, sep = "\t", header = F)
h_Cg   <- h_Cg[order(h_Cg$V2, decreasing = T),]
h_Cg_1 <- mapIds(x = org.Hs.eg.db, keys = h_Cg[1:i,1], keytype = "SYMBOL", column = "ENTREZID")
h_Cg_1 <- na.omit(h_Cg_1)
#开始整理数据：Healthy组别的基因
h_Hg   <- read.csv(file = "healthy_sortGene.rnk",stringsAsFactors = F, sep = "\t", header = F)
h_Hg   <- h_Hg[order(h_Hg$V2, decreasing = T),]
h_Hg_1 <- mapIds(x = org.Hs.eg.db, keys = h_Hg[1:i,1], keytype = "SYMBOL", column = "ENTREZID")
h_Hg_1 <- na.omit(h_Hg_1)
#开始整理数据：从UNIPROT随机抽取基因
set.seed(111)
h_Rg   <- base::sample(x = keys(org.Hs.eg.db, keytype = "UNIPROT"), size = i, replace = F)
h_Rg_1 <- mapIds(x = org.Hs.eg.db, keys = h_Rg, keytype = "UNIPROT", column = "ENTREZID")
h_Rg_1 <- na.omit(h_Rg_1)


#开始跑KEGG：智人的Cancer组别：跑KEGG----p值排序----取TOP
enrich.KEGG.h_Cg_1 = enrichKEGG(gene = h_Cg_1, organism = "hsa", keyType = "kegg",
                                pvalueCutoff = 0.05)
hC_kegg <- enrich.KEGG.h_Cg_1@result[order(enrich.KEGG.h_Cg_1@result$pvalue, decreasing = F),]
hC_kegg <- hC_kegg[1:10,]
#开始跑KEGG：智人的Healthy组别：跑KEGG----p值排序----取TOP
enrich.KEGG.h_Hg_1 = enrichKEGG(gene = h_Hg_1, organism = "hsa", keyType = "kegg",
                                pvalueCutoff = 0.05)
hH_kegg <- enrich.KEGG.h_Hg_1@result[order(enrich.KEGG.h_Hg_1@result$pvalue, decreasing = F),]
hH_kegg <- hH_kegg[1:10,]
#开始跑KEGG：智人的Random组别：跑KEGG----p值排序----取TOP
enrich.KEGG.h_Rg_1 = enrichKEGG(gene = h_Rg_1, organism = "hsa", keyType = "kegg",
                                pvalueCutoff = 0.05)
hR_kegg <- enrich.KEGG.h_Rg_1@result[order(enrich.KEGG.h_Rg_1@result$pvalue, decreasing = F),]
hR_kegg <- hR_kegg[1:10,]


#开始作图：智人的Cancer组别KEGG结果的TOP10
hC_kegg$Description <- factor(hC_kegg$Description,levels = rev(hC_kegg$Description))
p_1 <- ggplot(hC_kegg,aes(x=-log10(pvalue), y=Description))+
       geom_bar(stat = "identity",fill="black",width = 0.8) +
       theme_classic()+
       theme(text=element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.title = element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.text = element_text(size = 12,family = "serif",color = "black"))+
       ylab(NULL)
p_1; getwd()
ggsave(plot = p_1, filename = paste0(i,"erich.kegg.bar.cancer.hsa.jpg"), units = "cm", width = 14, height = 7, dpi = 800)

#开始作图：智人的Healthy组别KEGG结果的TOP10
hH_kegg$Description <- factor(hH_kegg$Description,levels = rev(hH_kegg$Description))
p_2 <- ggplot(hH_kegg,aes(x=-log10(pvalue), y=Description))+
       geom_bar(stat = "identity",fill="black",width = 0.8) +
       theme_classic()+
       theme(text=element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.title = element_text(size = 12,family = "serif",color = "black")) +
       theme(axis.text = element_text(size = 12,family = "serif",color = "black")) +
       ylab(NULL)
p_2; getwd()
ggsave(plot = p_2, filename = paste0(i,"erich.kegg.bar.healthy.hsa.jpg"), units = "cm", width = 14, height = 7, dpi = 800)

#开始作图：智人的Healthy组别KEGG结果的TOP10
hR_kegg$Description <- factor(hR_kegg$Description, levels = rev(hR_kegg$Description))
p_3 <- ggplot(hR_kegg,aes(x=-log10(pvalue), y=Description))+
       geom_bar(stat = "identity",fill="black",width = 0.8) +
       theme_classic()+
       theme(text=element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.title = element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.text = element_text(size = 12,family = "serif",color = "black"))+
       ylab(NULL)
p_3; getwd()
ggsave(plot = p_3, filename = paste0(i,"erich.kegg.bar.random.hsa.jpg"), units = "cm", width = 14, height = 7, dpi = 800)

#开始作图：合并Cancer、Healthy、Random组别
hU_kegg <- rbind(hC_kegg,hH_kegg,hR_kegg)
hU_kegg$ssot <- c(rep("cancer",10),rep("healthy",10),rep("random",10))
hU_kegg$Description <- factor(hU_kegg$Description, levels = rev(hU_kegg$Description))

p_4 <- ggplot(hU_kegg)+
       geom_bar(aes(x = -log10(pvalue), y = Description, fill = ssot), 
                stat = "identity", width = 0.8) +
       theme_classic()+
       theme(text=element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.title  = element_text(size = 12,family = "serif",color = "black"))+
       theme(axis.text   = element_text(size = 12,family = "serif",color = "black"),
             axis.text.y = element_text(angle = 0),
             plot.margin = unit(c(0.2,0.2,0.2,0.2), units = "cm"))+
       scale_fill_manual(breaks = hU_kegg$ssot,
                         values = c(rep("indianred3",10),rep("royalblue4",10),rep("gray40",10))) +
       theme(legend.position = c(0.8,0.1), 
             legend.text = element_text(size = 12,family = "serif",color = "black")) +
       guides(fill = guide_legend(title = "Group")) +
       #geom_vline(xintercept = -log10(0.05), linetype = 2, color = "red", size = 0.5) +
       ylab(NULL)
p_4; getwd()
ggsave(plot = p_4, filename = paste0(i,"erich.kegg.bar.UNION.hsa.jpg"), units = "cm", width = 18, height = 18, dpi = 800)
ggsave(plot = p_4, filename = paste0(i,"erich.kegg.bar.UNION.hsa.pdf"), units = "cm", width = 18, height = 18)















#安装包
#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("clusterProfiler")  #用来做富集分析
#BiocManager::install("topGO")  #画GO图用的
#BiocManager::install("Rgraphviz")
#BiocManager::install("pathview") #看KEGG pathway的
#BiocManager::install("org.Hs.eg.db") #这个包里存有人的注释文件


#AnnotationDbi::数据的面貌
#columns(org.Hs.eg.db)
#org.Hs.eg.db$GENENAME
#keys(org.Hs.eg.db, keytype = "UNIPROT")


#开始跑GO：举例
#erich.go.BP.down = enrichGO(gene = t_gene_up,
#                            OrgDb = org.Mm.eg.db,
#                            keyType = "ENTREZID",
#                            ont = "BP",
#                            pvalueCutoff = 0.05,
#                            qvalueCutoff = 0.05)
#C_kegg <- erich.go.BP.down@result[1:15,]



