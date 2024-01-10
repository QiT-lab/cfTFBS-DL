library(ChIPseeker)
library(ggplot2)
library(ggupset)
library(ggimage)
library(Guitar)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

cf_pred <- read.table('./result/prediction/FOXA1_breast_bc_02_predict.txt',sep = '\t',header = T)
cf_seq <- import_fst('./result/seq/bc_02_cf_seq.fst')
index <- which(cf_pred$score>=0.98)
bed <- cf_seq[index,]
bed$position <- gsub('>','',bed$position)
bed <- separate(bed,'position',into = c('chrom','start','end'),sep='[:-]')
write.table(bed[,1:3],'./result/FOXA1_breast_bc_02_peak.bed',sep = '\t',quote = F,row.names = F,col.names = F)

stBedFiles <- list("/mnt/data1/qiting/cfDNA/TF_binding/result/peak//CTCF_lung_G2_peak.bed",'/mnt/data1/qiting/cfDNA/TF_binding/result/peak/CTCF_bc_02_peak.bed','./result/peak/FOXA1_breast_bc_02_peak.bed')
p <- GuitarPlot(txTxdb = txdb, 
                stBedFiles = stBedFiles, 
                headOrtail = FALSE,
                enableCI = FALSE, 
                mapFilterTranscript = TRUE, 
                pltTxType = c("mrna"), 
                #stGroupName = c('CTCF_lung','CTCF_breast','FOXA1_breast')
                )

p <- p + ggtitle(label = "TFBS Distribution on mRNA") + 
  xlab("") +  
  ylab("Mean Depth") +
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(panel.background = element_blank(),
                    panel.grid.major = element_blank(),
                    axis.line = element_line(size = 0.5, colour = 'black'), 
                    axis.title = element_text(size=13),
                    axis.text = element_text(size = 10),
                    plot.title = element_text(size = 16,hjust = 0.5),
                    legend.position = c(0.3,0.9),
                    legend.background = element_blank(),
                    legend.key = element_blank(),
                    legend.key.size = unit(5,'mm'),
                    legend.title = element_blank())
p
ggsave('./result/peak/lung_bc_distribution_on_mRNA.pdf',width = 6,height = 4)

f <- readPeakFile(stBedFiles[[1]])
#f <- GRanges(seqnames = bed$chrom,IRanges(start = bed$start,end = bed$end))
covplot(f)
ggsave('./result/peak/FOXA1_bc_02_peaks_chromosomes.pdf',width = 10,height = 10)

x = annotatePeak(f, tssRegion=c(-3000, 3000), TxDb=txdb,ignoreUpstream = T,ignoreDownstream = T)
plotAnnoPie(x)
ggsave('./result/peak/CTCF_lung_G2_peak_pie.pdf',width=6,height = 4)
plotAnnoBar(x)
plotAnnoPie.csAnno(x)
plotDistToTSS(x)
ggsave('./result/peak/FOXA1_bc_02_peak_TSS_bar.pdf',width=6,height = 4)
upsetplot(x, vennpie=TRUE)


peakHeatmap(f,TxDb=txdb, upstream=3000, downstream=3000)


require(clusterProfiler)
# 将bed文件读入（readPeakFile是利用read.delim读取，然后转为GRanges对象）
seq=lapply(stBedFiles[[1]], readPeakFile)

genes=lapply(seq, function(i) seq2gene(i, c(-1000, 3000), 3000, TxDb=txdb))
go <- enrichGO(gene=genes[[1]],OrgDb='org.Hs.eg.db',ont= "CC",pAdjustMethod="BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
dotplot(go)
ggsave('./result/peak/FOXA1_bc_02_peak_gene_go.pdf',width=6,height = 4)
