library(tidyverse)
library(data.table)
library(tidyfst)
#install_github("remap-cisreg/ReMapEnrich")
library(ReMapEnrich)
library(Biostrings)
library(keras)
library(mltools)
library(ggseqlogo)
library(tidyfst)
library(pROC)
library(ChIPseeker)
library(ggplot2)
library(ggupset)
library(ggimage)
library(Guitar)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

my_theme <- theme(panel.background = element_blank(),
                  panel.grid.major = element_blank(),
                  axis.line = element_line(size = 0.5, colour = 'black'), 
                  axis.title = element_text(size=13),
                  axis.text = element_text(size = 10),
                  plot.title = element_text(size = 16,hjust = 0.5),
                  legend.position = c(0.9,0.9),
                  legend.background = element_blank(),
                  legend.key = element_blank(),
                  legend.key.size = unit(5,'mm'),
                  legend.title = element_blank())


oneHotEncode <- function(bio_seq)
{
  test <- data.frame(row.names =c('A','T','C','G'), col1=c(1,0,0,0),col2=c(0,1,0,0),col3=c(0,0,1,0),col4=c(0,0,0,1))
  bio_char <- strsplit(bio_seq,'')[[1]]
  df <- unlist(lapply(1:length(bio_char), function(i) {return(test[bio_char[i],])}))
  return(df)
}

getfasta <- function(bed)
{
  write.table(bed,paste0(dir,'test_peak.bed'),sep = '\t',quote = F,row.names = F,col.names = F)
  cmd <- paste('/home/qiting/software/miniconda/miniconda3/bin/bedtools getfasta -fi /mnt/public/reference/human/hg38/bwa/hg38.fa -bed',paste0(dir,'test_peak.bed'),'-fo',paste0(dir,'peak_seq.fa'),sep = ' ')
  system(cmd)
  seq <- fread(paste0(dir,'peak_seq.fa'),header = F)
  s <- seq[seq(1,nrow(seq),2),]
  d <- seq[seq(2,nrow(seq),2),]
  seq <- cbind(s,d)
  names(seq) <- c('position','sequence')
  seq$sequence <- toupper(seq$sequence)
  return(seq)
}

calculate_gc <- function(sequence) {
  gc_content <- sum(strsplit(sequence, "")[[1]] %in% c("G", "C")) / nchar(sequence)
  return(gc_content)
}

get_samelength <- function(aa)
{
  ll <- length(aa)
  if(ll < 276){
    aa <- c(aa,rep(0,(276-ll)))
  }
  return(aa)
}
