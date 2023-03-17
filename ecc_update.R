#Write by Sun Haiyang on 2023.02.27
#01、the Environment
{
  rm(list=ls())
  for (i in c("data.table","Repitools","stringr","Biostrings","tidyr","liftOver","enrichplot","patchwork","ggpubr","ggsci","ggtranscript","rtracklayer","BiocManager","dplyr","IRanges","GenomicRanges","ReactomePA","ggplot2","clusterProfiler","topGO","Rgraphviz","devtools","BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Mmusculus.UCSC.mm9","BSgenome.Mmusculus.UCSC.mm10","org.Mm.eg.db")){library(i, character.only = T,quietly = T)}; rm(i)
}
#02、yield Index
{
  path = "/disk5/shy/1shyJOB/7eccdna/aaUpdate/"
  ##yield BED of specifical type of genes
  {
    shy_yield_allGen <- function(inpuFile, saveFile){
      thePath="/Reference/aaSHT/eccDNA/1GTFGFF/hg38/"
      posi <- import(paste0(thePath,inpuFile))
      posi <- as.data.frame(posi[which(posi$type=="gene")])[,c("seqnames","start","end","gene_name")]
      posi <- posi[!duplicated(posi$gene_name),]
      write.table(file = paste0(thePath,saveFile),
                  x = posi,sep = "\t", quote = F, col.names = F, row.names = F)
    }
    shy_yield_allGen(inpuFile = "gencode.v21.annotation.gtf", saveFile = "allGene.hg38.position")
    shy_yield_proten <- function(inpuFile, saveFile){
      thePath="/Reference/aaSHT/eccDNA/1GTFGFF/hg38/"
      posi <- import(paste0(thePath,inpuFile))
      posi <- as.data.frame(posi[which(posi$type=="gene" & posi$gene_type=="protein_coding")])[,c("seqnames","start","end","gene_name")]
      posi <- posi[!duplicated(posi$gene_name),]
      write.table(file = paste0(thePath,saveFile),
                  x = posi,sep = "\t", quote = F, col.names = F, row.names = F)
    }
    shy_yield_proten(inpuFile = "gencode.v21.annotation.gtf", saveFile = "protein_coding.hg38.position")
    shy_yield_ncRNAs <- function(inpuFile, saveFile){
      thePath="/Reference/aaSHT/eccDNA/1GTFGFF/hg38/"
      posi <- import(paste0(thePath,"gencode.v21.annotation.gtf"))
      posi <- as.data.frame(posi[which(posi$type=="gene")])
      posi <- posi[grep(posi$gene_type,pattern = "*RNA*"),]
      posi <- posi[which(posi$gene_type!="lincRNA"),c("seqnames","start","end","gene_name","gene_type")]
      posi <- posi[!duplicated(posi$gene_name),]
      posi$thename <- paste0(posi$gene_name,"(",posi$gene_type,")")
      write.table(file = paste0(thePath,saveFile),
                  x = posi[,c(1,2,3,6)],sep = "\t", quote = F, col.names = F, row.names = F)
    }
    shy_yield_ncRNAs(inpuFile = "gencode.v21.annotation.gtf", saveFile = "ncRNA.hg38.position")
    shy_yield_lncRNA <- function(inpuFile, saveFile){
      thePath="/Reference/aaSHT/eccDNA/1GTFGFF/hg38/"
      posi <- import(paste0(thePath,inpuFile))
      posi <- as.data.frame(posi[which(posi$type=="gene")])[,c("seqnames","start","end","gene_name")]
      posi <- posi[!duplicated(posi$gene_name),]
      write.table(file = paste0(thePath,saveFile),
                  x = posi,sep = "\t", quote = F, col.names = F, row.names = F)
    }
    shy_yield_lncRNA(inpuFile = "gencode.v21.long_noncoding_RNAs.gtf", saveFile = "lncRNA.hg38.position")
    shy_yield_pseudo <- function(inpuFile, saveFile){
      thePath="/Reference/aaSHT/eccDNA/1GTFGFF/hg38/"
      posi <- import(paste0(thePath,inpuFile))
      posi <- as.data.frame(posi)[,c("seqnames","start","end","gene_id")]
      posi <- posi[!duplicated(posi$gene_id),]
      write.table(file = paste0(thePath,saveFile),
                  x = posi,sep = "\t", quote = F, col.names = F, row.names = F)
    }
    shy_yield_pseudo(inpuFile = "gencode.v21.2wayconspseudos.gtf", saveFile = "pseudos.hg38.position")
  }
}
#03：data of Version 1
##based on hg38 + length shoter than 30000000
{
shy_obtain_position <- function(release="hg38", cutoff=30000000){
  tm_1 <- Sys.time()
  ecc1 <- fread("/disk5/shy/1shyJOB/7eccdna/raw/hsa.txt", sep="\t")
  raws <- ecc1[which(ecc1$Length < cutoff),c(1,3,5,6,7,48,49,50,51)]
  raws <- raws[!(is.na(raws$`Coordinate(hg38)`)),]
  raws <- raws[grep(raws$`Coordinate(hg38)`, pattern = "*KI|*GL|_", invert = T), ]
  print(paste0("fasta file size is about ",round(sum(raws$Length)/(1024*1024*1024),1)," Gb"))
  pov1 <- tidyr::separate(raws[, c(1,2)], col = paste0("Coordinate(",release,")"), into = c("seqnames","start","end"),sep = ":|-",remove = T)
  raws <<- raws
  pov1 <<- pov1
  write.table(pov1[,c(2,3,4,1)], file = paste0(path,"bedV1/",release,".bed"), col.names = F, quote = F, sep = "\t", row.names = F)
  tm_2 <- Sys.time()
  print(tm_2-tm_1)
}
shy_obtain_position()
}
{
  shell <- paste0(path, "src/make_beds.sh")
  cat("#!/bin/bash\n", file = shell)
  bedf <- paste0(path,"bedV1/",release,".bed")
  cmd_31 <- paste0("sort-bed ",bedf," >",paste0(path,"bedV1/",release,".sort.bed\n"))
  cmd_32 <- paste0("cat ",bedf," |head -n 580000 >", paste0(path,"bedV1/",release,".bd_1.bed\n"))
  x=580000; y=20000
  options(scipen = 200)
  for (i in seq(round((720000-580000)/20000, digits = 1))){
    x=x+y
    assign(paste0("cmd_3",i+2),paste0("cat ",bedf," |head -n ",x," |tail -n ",y, " >", paste0(path,"bedV1/",release,".bd_",i+1,".bed\n")))
    }
  cmd_3a <- paste0("cat ",bedf," |tail -n 9369 >", paste0(path,"bedV1/",release,".bd_9.bed\n"))
  for (i in cmd_31) {cat(i, file = shell, append = T)}
  for (i in cmd_32) {cat(i, file = shell, append = T)}
  for (i in cmd_33) {cat(i, file = shell, append = T)}
  for (i in cmd_34) {cat(i, file = shell, append = T)}
  for (i in cmd_35) {cat(i, file = shell, append = T)}
  for (i in cmd_36) {cat(i, file = shell, append = T)}
  for (i in cmd_37) {cat(i, file = shell, append = T)}
  for (i in cmd_38) {cat(i, file = shell, append = T)}
  for (i in cmd_39) {cat(i, file = shell, append = T)}
  for (i in cmd_3a) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"src_log/"),basename(shell),".log", " 2>&1 &"))
}
#04：statistics
#unique(unname(unlist(sapply(unique(raws$Disease), function(x){str_split(x, pattern = ", ")}))))
#unique(unname(unlist(sapply(unique(raws$`Cell/Tissues`), function(x){str_split(x, pattern = ", ")}))))
#unique(unname(unlist(sapply(unique(raws$Identify_method), function(x){str_split(x, pattern = ", ")}))))
#
#
#
#
#05：The Annotation Pipeline
{
  shell <- paste0(path,"src/annotationPipeline.sh")
  cat("#!/bin/bash\n", file = shell)
  bed1 <- list.files(path = "/disk5/shy/1shyJOB/7eccdna/aaUpdate/bedV1", pattern = "hg38.bed", full.names = T)
  bed2 <- list.files(path = "/disk5/shy/1shyJOB/7eccdna/aaUpdate/bedV1", pattern = "*sort*", full.names = T)
  bed3 <- list.files(path = "/disk5/shy/1shyJOB/7eccdna/aaUpdate/bedV1", pattern = "*bd_*", full.names = T)
  an00 <- list.files(path = "/Reference/aaSHT/appTools/liftOver", pattern = "hg38*ToHg19*", full.names = T)
  an01 <- list.files(path = "/Reference/aaSHT/eccDNA/1GTFGFF/hg38", pattern = "*position", full.names = T)
  an02 <- list.files(path = "/Reference/aaSHT/eccDNA/3repeat", pattern = "hg38*", full.names = T)
  an03 <- list.files(path = "/Reference/aaSHT/eccDNA/4tss", pattern = "hg38*", full.names = T)
  an04 <- list.files(path = "/Reference/aaSHT/eccDNA/5utr", pattern = "hg38*", full.names = T)
  an05 <- list.files(path = "/Reference/aaSHT/eccDNA/6exon", pattern = "hg38*", full.names = T)
  an06 <- list.files(path = "/Reference/aaSHT/eccDNA/7intron", pattern = "hg38*", full.names = T)
  an07 <- list.files(path = "/Reference/aaSHT/eccDNA/14seq", pattern = "hg38.fa$", full.names = T)
  an08 <- list.files(path = "/Reference/aaSHT/eccDNA/8cpg", pattern = "hg38*", full.names = T)
  an09 <- list.files(path = "/Reference/aaSHT/eccDNA/9promoter", pattern = "hg38*", full.names = T)
  an10 <- list.files(path = "/Reference/aaSHT/eccDNA/10enhancer", pattern = "hg38*", full.names = T)
  an11 <- list.files(path = "/Reference/aaSHT/eccDNA/11tfbs", pattern = "hg38*", full.names = T)
  an12 <- list.files(path = "/Reference/aaSHT/eccDNA/12histone", pattern = "hg38*", full.names = T)
  an13 <- list.files(path = "/Reference/aaSHT/eccDNA/13DNaseI", pattern = "hg38*", full.names = T)
  ##hg19
  cmd_a1 <- paste0("liftOver ", bed1," ",an00," ",paste0(path,"src_Annotation/hg19.liftOver.bed "),paste0(path,"src_Annotation/hg19.unmapped.bed\n"))
  cmd_a2 <- paste0("mv ",paste0(path,"src_Annotation/hg19.liftOver.bed "),paste0(path,"src_Annotation/a1_hg19.liftOver.txt\n"),
                   "rm ",paste0(path,"src_Annotation/hg19.unmapped.bed\n"))
  ##genes
  cmd_a3 <- paste0("bedtools intersect -a ", bed1," -b ",an01," -wa -wb -loj |","sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,collapse |",
                   "sed 's/1\t\\./0\t\\./g; s/,/, /g' |cut -f 4,5,6 >",
                   paste0(path,"src_Annotation/","a",seq(length(an01)+1, from=2),"_",str_replace(basename(an01),"position","bedtools.txt\n")))
  cmd_a4 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed1," ",an01,
                   " |cut -f 4,5 >",paste0(path,"src_Annotation/","a",seq(length(an01)+1, from=2),"_",str_replace(basename(an01),"position","bedmap.txt\n")))
  ##repeats
  cmd_a5 <- paste0("bedtools intersect -a ", bed3," -b ",an02," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g; s/,/, /g' |cut -f 4,5,6 >",
                   paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  cmd_a6 <- paste0("cat ",paste(paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt")), collapse = " "),
                   " >",paste0(path,"src_Annotation/","a7_",str_replace(basename(an02),"bed","bedtools.txt\n")))
  cmd_a7 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  cmd_a8 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed2," ",an02,
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/","a7_",str_replace(basename(an02),"bed","bedmap_pre.txt\n")))
  cmd_a9 <- paste0("sort -k 1 ",
                   paste0(path,"src_Annotation/","a7_",str_replace(basename(an02),"bed","bedmap_pre.txt >")),
                   paste0(path,"src_Annotation/","a7_",str_replace(basename(an02),"bed","bedmap.txt\n")))
  cmd_ab <- paste0("rm ",paste0(path,"src_Annotation/","a7_",str_replace(basename(an02),"bed","bedmap_pre.txt\n")))
  ##TSS
  cmd_b1 <- paste0("bedtools intersect -a ", bed1," -b ",an03," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g' |cut -f 4,5 >",
                   paste0(path,"src_Annotation/","b1_",str_replace(basename(an03),"bed","bedtools.txt\n")))
  ##UTRs
  cmd_b2 <- paste0("bedtools intersect -a ", bed1," -b ",an04," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g' |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",c("b3_","b2_"),str_replace(basename(an04),"bed","bedtools.txt\n")))
  cmd_b3 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed1," ",an04,
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",c("b3_","b2_"),str_replace(basename(an04),"bed","bedmap.txt\n")))
  ##Exon & Intron
  cmd_b4 <- paste0("bedtools intersect -a ", bed1," -b ",c(an05,an06)," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g' |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",c("b4_","b5_"),str_replace(basename(c(an05,an06)),"bed","bedtools.txt\n")))
  cmd_b5 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed2," ",c(an05,an06),
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",c("b4_","b5_"),str_replace(basename(c(an05,an06)),"bed","bedmap_pre.txt\n")))
  cmd_b6 <- paste0("sort -k 1 ",
                   paste0(path,"src_Annotation/",c("b4_","b5_"),str_replace(basename(c(an05,an06)),"bed","bedmap_pre.txt >")),
                   paste0(path,"src_Annotation/",c("b4_","b5_"),str_replace(basename(c(an05,an06)),"bed","bedmap.txt\n")))
  cmd_b7 <- paste0("rm ",paste0(path,"src_Annotation/",c("b4_","b5_"),str_replace(basename(c(an05,an06)),"bed","bedmap_pre.txt\n")))
  ##GC content
  cmd_c1 <- paste0("bedtools nuc -fi ", an07," -bed ",bed3," |grep -v gc |cut -f 4,6 >",
                   paste0(path,"src_Annotation/","c1_",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  cmd_c2 <- paste0("cat ",paste(paste0(path,"src_Annotation/","c1_",str_replace(basename(bed3),"bed","bedtools.txt")), collapse = " "),
                   " >",paste0(path,"src_Annotation/","c1_",str_replace(basename(an07),"fa","gc.bedtools.txt\n")))
  cmd_c3 <- paste0("rm ",paste0(path,"src_Annotation/","c1_",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  ##CpG island
  cmd_c4 <- paste0("bedtools intersect -a ", bed1," -b ",an08," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g' |cut -f 4,5 >",
                   paste0(path,"src_Annotation/","c2_",str_replace(basename(an08),"bed","bedtools.txt\n")))
  cmd_c5 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed2," ",an08,
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",str_replace(basename(an08),"bed","bedmap_pre.txt\n")))
  cmd_c6 <- paste0("sort -k 1 ",
                   paste0(path,"src_Annotation/",str_replace(basename(an08),"bed","bedmap_pre.txt >")),
                   paste0(path,"src_Annotation/","c2_",str_replace(basename(an08),"bed","bedmap.txt\n")))
  cmd_c7 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(an08),"bed","bedmap_pre.txt\n")))
  ##Promoter
  cmd_d1 <- paste0("bedtools intersect -a ", bed1," -b ",an09," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g' |cut -f 4,5,6 >",
                   paste0(path,"src_Annotation/","d1_",str_replace(basename(an09),"bed","bedtools.txt\n")))
  cmd_d2 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed2," ",an09,
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",str_replace(basename(an09),"bed","bedmap_pre.txt\n")))
  cmd_d3 <- paste0("sort -k 1 ",
                   paste0(path,"src_Annotation/",str_replace(basename(an09),"bed","bedmap_pre.txt >")),
                   paste0(path,"src_Annotation/","d1_",str_replace(basename(an09),"bed","bedmap.txt\n")))
  cmd_d4 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(an09),"bed","bedmap_pre.txt\n")))
  ##Enhancer
  cmd_d5 <- paste0("bedtools intersect -a ", bed1," -b ",an10," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g' |cut -f 4,5 >",
                   paste0(path,"src_Annotation/","d2_",str_replace(basename(an10),"bed","bedtools.txt\n")))
  cmd_d6 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed2," ",an10,
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",str_replace(basename(an10),"bed","bedmap_pre.txt\n")))
  cmd_d7 <- paste0("sort -k 1 ",
                   paste0(path,"src_Annotation/",str_replace(basename(an10),"bed","bedmap_pre.txt >")),
                   paste0(path,"src_Annotation/","d2_",str_replace(basename(an10),"bed","bedmap.txt\n")))
  cmd_d8 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(an10),"bed","bedmap_pre.txt\n")))
  ##TFBS
  cmd_e1 <- paste0("bedtools intersect -a ", bed3," -b ",an11," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g; s/,/, /g' |cut -f 4,5,6 >",
                   paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  cmd_e2 <- paste0("cat ",paste(paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt")), collapse = " "),
                   " >",paste0(path,"src_Annotation/","e1_",str_replace(basename(an11),"bed","bedtools.txt\n")))
  cmd_e3 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  ##Histone
  cmd_e4 <- paste0("bedtools intersect -a ", bed3," -b ",an12," -wa -wb -loj |",
                   "sed 's/\t\t/\t/g' |",
                   "bedtools groupby -i - -g 1-4 -c 8 -o count,distinct |",
                   "sed 's/1\t\\./0\t\\./g; s/,/, /g' |cut -f 4,5,6 >",
                   paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  cmd_e5 <- paste0("cat ",paste(paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt")), collapse = " "),
                   " >",paste0(path,"src_Annotation/","e2_",str_replace(basename(an12),"bed","bedtools.txt\n")))
  cmd_e6 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(bed3),"bed","bedtools.txt\n")))
  ##DNaseI
  cmd_e7 <- paste0("bedmap --faster --echo --bases-uniq-f --delim ","'\t' ", bed2," ",an13,
                   " |cut -f 4,5 >",
                   paste0(path,"src_Annotation/",str_replace(basename(an13),"bed","bedmap_pre.txt\n")))
  cmd_e8 <- paste0("sort -k 1 ",
                   paste0(path,"src_Annotation/",str_replace(basename(an13),"bed","bedmap_pre.txt >")),
                   paste0(path,"src_Annotation/","e3_",str_replace(basename(an13),"bed","bedmap.txt\n")))
  cmd_e9 <- paste0("rm ",paste0(path,"src_Annotation/",str_replace(basename(an13),"bed","bedmap_pre.txt\n")))
  for (i in c(cmd_a1,cmd_a2,cmd_a3,cmd_a4,cmd_a5,cmd_a6,cmd_a7,cmd_a8,cmd_a9,
              cmd_ab,cmd_b1,cmd_b2,cmd_b3,cmd_b4,cmd_b5,cmd_b6,cmd_b7,cmd_c1,
              cmd_c2,cmd_c3,cmd_c4,cmd_c5,cmd_c6,cmd_c7,cmd_d1,cmd_d2,cmd_d3,
              cmd_d4,cmd_d5,cmd_d6,cmd_d7,cmd_e1,cmd_e2,cmd_e3,cmd_e4,cmd_e5,
              cmd_e6,cmd_e7,cmd_e8,cmd_e9)) {cat(i, file = shell, append = T)}
  print(paste0("nohup bash ",shell, " >",paste0(path,"src_log/"),basename(shell),".log", " 2>&1 &"))
  rm(cmd_a1,cmd_a2,cmd_a3,cmd_a4,cmd_a5,cmd_a6,cmd_a7,cmd_a8,cmd_a9,cmd_ab,
     cmd_b1,cmd_b2,cmd_b3,cmd_b4,cmd_b5,cmd_b6,cmd_b7,
     cmd_c1,cmd_c2,cmd_c3,cmd_c4,cmd_c5,cmd_c6,cmd_c7,
     cmd_d1,cmd_d2,cmd_d3,cmd_d4,cmd_d5,cmd_d6,cmd_d7,cmd_d8,
     cmd_e1,cmd_e2,cmd_e3,cmd_e4,cmd_e5,cmd_e6,cmd_e7,
     bed1,bed2,bed3,
     an01,an02,an03,an04,an05,an06,an07,an08,an09,an10,an11,an12,an13,an00)
}
{
  leng <- as.numeric(pov1$end)-as.numeric(pov1$start)
  write.table(x = data.frame(idx = pov1$eccDNA_ID, length = leng),
              col.names = F, row.names = F, quote = F,sep = "\t",
              file = paste0(path,"src_Annotation/", "f1_hg38.Length.txt"))
  hg38 <- paste0(pov1$seqnames,":",pov1$start,"-", pov1$end)
  write.table(x = data.frame(idx = pov1$eccDNA_ID, length = hg38),
              col.names = F, row.names = F, quote = F, sep = "\t",
              file = paste0(path,"src_Annotation/", "f1_hg38.Coordinate.txt"))
}
#
#
#
#
#06：to Process files
{
  file <- read.table(file = paste0(path,"src_Annotation/a1_hg19.liftOver.txt"), header = F, sep = "\t", stringsAsFactors = F)
  file$co19 <- paste0(file$V1,":",file$V2,"-",file$V3)
  file <- file[,c(4,5)]
  tm19 <- rbind(file, data.frame(V4 = dplyr::setdiff(pov1$eccDNA_ID, file$V4),
                                 co19 = rep("NA", length(dplyr::setdiff(pov1$eccDNA_ID, file$V4)))))
  tm19 <- tm19[order(tm19$V4, decreasing = F),]
  hg19 <- tm19$co19
  rm(tm19, file)
}



