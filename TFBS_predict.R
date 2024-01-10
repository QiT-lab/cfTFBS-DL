library(tidyverse)
library(data.table)
#install_github("remap-cisreg/ReMapEnrich")
library(ReMapEnrich)
library(Biostrings)
library(keras)


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


gen_neg_seq <- function(chromosomes,chrom_size,temp,width,seq)
{
  gc <- temp$gc
  random_gc=0
  neg_seq <- 'NNNN'
  while(length(grep('N',neg_seq))!=0 | abs(random_gc-gc) > 0.1 ){
    random_chrom <- sample(chromosomes,1)
    chrom_length <- chrom_size[which(chrom_size$V1==random_chrom),2]
    random_start <- sample(1:chrom_length-width,1)
    random_end <- random_start+width
    same_chr <- seq[which(seq$chrom==random_chrom),]
    t <- apply(same_chr, 1, function(x) {intersect(random_start:random_end,x[2]:x[3])})
    if(length(t)==0){ 
      cmd <- paste('samtools faidx /mnt/public/reference/human/hg38/bwa/hg38.fa',paste0(random_chrom,':',random_start,'-',random_end),sep = ' ')
      neg_seq <- system(cmd,intern = TRUE)
      position <- neg_seq[1]
      neg_seq <- toupper(paste0(neg_seq[2],neg_seq[3]))
      random_gc <- calculate_gc(neg_seq)
    }
  }
  temp$neg_position <- position
  temp$neg_sequence <- neg_seq
  temp$neg_gc <- random_gc
  return(temp)
}





# peak_g <- bedToGranges('/mnt/data1/qiting/cfDNA/TF_binding/data/remap2022_peripheral-blood-mononuclear-cell_all_macs2_hg38_v1_0.bed')
# shuffle_seq <- shuffle(peak_g,byChrom = TRUE)
# shuf <- resize(shuffle_seq,51,fix = 'center')
# shuf <- as.data.frame(shuf)
# shuf$start <- shuf$start-1
# shuf_seq <- getfasta(shuf[,1:3])
# shuf_seq <- shuf_seq[-grep('N',shuf_seq$sequence),]
# shuf_seq$label <- '0'
# 
# set.seed(11)
# seq <- seq[sample(1:nrow(seq),nrow(shuf_seq)),]
# seq$label <- '1'


### one-hot coding
#library(hashmap,lib.loc = '/home/qiting/software/miniconda/miniconda3/envs/hashmap/lib/R/library/')
oneHotEncode <- function(bio_seq)
{
  test <- data.frame(row.names =c('A','T','C','G'), col1=c(1,0,0,0),col2=c(0,1,0,0),col3=c(0,0,1,0),col4=c(0,0,0,1))
  bio_char <- strsplit(bio_seq,'')[[1]]
  df <- unlist(lapply(1:length(bio_char), function(i) {return(test[bio_char[i],])}))
  return(df)
}

cnn_model <- function(shuffle_seq,pos_matrix,neg_matrix,q)
{
  set.seed(12)
  test <- sample(1:nrow(shuffle_seq),q,replace = F)
  train <- setdiff(1:nrow(shuffle_seq),test)
  
  test_data <- rbind(pos_matrix[test,],neg_matrix[test,])
  test_label <- c(rep(1,q),rep(0,q))
  
  train_data <- rbind(pos_matrix[train,],neg_matrix[train,])
  train_label <- c(rep(1,length(train)),rep(0,length(train)))
  
  s <- sample(1:nrow(train_data))
  train_data <- train_data[s,]
  train_label <- train_label[s]
  d <- sample(1:nrow(test_data))
  test_data <- test_data[d,]
  test_label <- test_label[d]
  
  model <- keras_model_sequential() %>%
    layer_conv_1d(filters =128, kernel_size = 76, activation = "relu", input_shape = c(324,1),strides = 4) %>%
    layer_max_pooling_1d(pool_size = 4) %>%
    layer_dropout(0.5) %>%
    layer_lstm(c(10),return_sequences = F) %>%
    #layer_flatten() %>%
    # layer_dense(units = 64, activation = "relu") %>%
    # #layer_lstm(c(20),return_sequences = F) %>%
    # layer_dense(units = 32, activation = "relu") %>%
    #layer_dense(units = 16, activation = "relu") %>%
    layer_dense(units = 1,activation = 'sigmoid')
  
  model %>% compile(optimizer = 'Adam',
                    loss = 'binary_crossentropy',
                    metrics = c('AUC'))
  
  history <- model %>% fit(train_data,train_label,epochs = 10,batch_size = 32,validation_split = 0.1)
  result <- model %>% evaluate(test_data,test_label)
  print(result)
  save_model_hdf5(model,'./model/CTCF_breast_model.hdf5')
  return(model)
}

###########################################################对chip-seq数据的分析
setwd('/mnt/data1/qiting/cfDNA/TF_binding/')
dir <- '/mnt/data1/qiting/cfDNA/TF_binding/data/CTCF_data/'
peak <- read.table('/mnt/data1/qiting/cfDNA/TF_binding/data/CTCF_data/ENCSR304XUZ.CTCF.breast_epithelium.bed')
peak <- peak[order(peak$V5,decreasing = T)[1:10000],]
names(peak) <- c('chrom','start','end','name','score','strand','peakStart','peakEnd','itemRgb')

bed <- unique(peak[,c(1,7,8)])
bed$peakStart <- bed$peakStart-40
bed$peakEnd <- bed$peakEnd+40
seq <- getfasta(bed)

seq$gc <- unlist(lapply(1:nrow(seq),function(i) {calculate_gc(seq$sequence[i])}))
seq$position <- gsub('>','',seq$position)
seq <- separate(seq,'position',into = c('chrom','start','end'),sep='[:-]')
seq$start <- as.integer(seq$start)
seq$end <- as.integer(seq$end)

chromosomes <- c(paste0("chr", 1:22),'chrX','chrY')
chrom_size <- read.table('./data/chromsize.txt')
chrom_size <- chrom_size[which(chrom_size$V1 %in% chromosomes),]
width=80

library(parallel)

detectCores()  
clus <- makeCluster(60)
clusterExport(clus,deparse(substitute(gen_neg_seq)))
clusterExport(clus,"chromosomes",envir = environment())
clusterExport(clus,"chrom_size",envir = environment())
clusterExport(clus,"seq",envir = environment())
clusterExport(clus,"width",envir = environment())
clusterExport(clus,deparse(substitute(calculate_gc)))

shuffle_seq <- parLapply(clus,1:nrow(seq), function(i) gen_neg_seq(chromosomes,chrom_size,seq[i,],width,seq))
shuffle_seq <- do.call(rbind,shuffle_seq)
#shuffle_seq <- shuffle_seq[-3722,]
stopCluster(clus)
write.table(shuffle_seq,'./data/CTCF_data/ENCSR304XUZ.CTCF.breast_shuffle_81.txt',quote = F,row.names = F)

shuffle_seq <- read.table('./data/GSE43098_CTCF_shuffle.txt',sep = '',header = T)
shuffle_seq <- shuffle_seq[-3722,]

pos_matrix <- do.call(rbind,lapply(shuffle_seq$sequence, oneHotEncode))
neg_matrix <- do.call(rbind,lapply(shuffle_seq$neg_sequence, oneHotEncode))

model <- cnn_model(shuffle_seq,pos_matrix,neg_matrix,q=1000)

pred <- model %>% predict(pos_matrix)
t <- which(pred>=1)

#############可视化卷积层的结果，得到seqlogo
library(ggseqlogo)
seq_matrix <- pos_matrix[t,]
sequence <- shuffle_seq$sequence[t]

motif_value <- get_seqlogo(seq_matrix,sequence,model,p=0.9)
win_value <- motif_value[[2]]
motif <- motif_value[[1]]
len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))
filter_index <- which(len_m>20)
ggseqlogo(motif[filter_index])

get_seqlogo <- function(seq_matrix, sequence,model,p)
{
  layer_outputs <- lapply(model$layers, function(layer) layer$output)
  activation_model <- keras_model(inputs = model$inputs,outputs = layer_outputs)
  activations <- activation_model %>% predict(seq_matrix)
  first_layer <- activations[[1]]
  tt <- list()
  win_value <- list()
  pos_len <- list()
  for (k in 1:128) {
    m <- first_layer[,,k]
    aa <- which(m>p,arr.ind = T)
    win_value[[k]] <- m[aa]
    pos_len[[k]] <- apply(aa,1,function(x) {nchar(sequence[x[1]])})
    tt[[k]] <- apply(aa,1,function(x) {substr(sequence[x[1]],x[2],x[2]+18)})
  }
  return(list(tt,win_value,pos_len))
}

##########################计算ROC
library(pROC)
library(ggplot2)
## roc的计算，可以一次性批量计算a、b、c三组数据
res1 <- roc(test_label,as.vector(pred_test),aur=TRUE,ci=TRUE)
res2 <- roc(train_label,as.vector(pred_train),aur=TRUE,ci=TRUE)
res3 <- roc(labels,as.vector(pred_all),aur=TRUE,ci=TRUE)
res <- list(`test: AUC=0.819`=res1,`train: AUC=0.875`=res2)
ggroc(res,legacy.axes = T)+
  #geom_abline(intercept = 0,slope=270)+
  annotate(geom = 'segment',x=0,y=0,xend = 1,yend = 1,lty=2)+
  #annotate(geom = 'text',label='AUC=0.8584',size=5,color='red',x=0.3,y=0.5)+
  my_theme+
  theme(panel.border = element_rect(size=1,fill = 'transparent'),
        legend.position = c(0.8,0.3))

##################################################################################
##find TFBS in cfDNA sequences

cfdna <- data.table::fread('/mnt/data1/qiting/cfDNA/TF_binding/data/cancer_data/bed/GSE171434_etana_h6_alexis_bc_02_sorted_center.bed',sep = '\t')
cf_bed <- unique(cfdna[,c(7:10,12)])
cf_bed <- cf_bed[which((cf_bed$V10<82) & (cf_bed$V10>30) & (cf_bed$V7 %in% chromosomes)),]
cf_bed <- cf_bed[,c(1,2,3)]
#cf_bed <- unique(cf_bed)
cf_seq <- getfasta(cf_bed)
#cf_seq$len <- nchar(cf_seq$sequence)
#cf_seq_test <- cf_seq[which(cf_seq$len==101),]


get_samelength <- function(aa)
{
  ll <- length(aa)
  if(ll < 324){
    aa <- c(aa,rep(0,324-ll))
  }
  return(aa)
}


get_cf_matrix <- function(temp_seq)
{
  library(parallel)
  detectCores()  
  clus <- makeCluster(60)
  clusterExport(clus,deparse(substitute(oneHotEncode)))
  #clusterExport(clus,"cf_seq",envir = environment())
  cf_list <- parLapply(clus,temp_seq$sequence, oneHotEncode)
  stopCluster(clus)
  
  cf_matrix <- do.call(rbind,lapply(1:length(cf_list), function(x) {get_samelength(cf_list[[x]])}))
  cf_matrix[which(is.na(cf_matrix),arr.ind = T)] <- 0
  return(cf_matrix)
}
aa <- c(seq(1,nrow(cf_seq),1000000),nrow(cf_seq)+1)
cf_matrix <- matrix(nrow = 0, ncol = 324)
for (i in 3:11) {
  print(i)
  temp_seq <- cf_seq[aa[i]:(aa[i+1]-1),]
  temp_matrix <- get_cf_matrix(temp_seq)
  cf_matrix <- rbind(cf_matrix,temp_matrix)
}
write.table(cf_matrix,'./data/cancer_data/bc_02_matrix.txt',quote = F,sep = '\t')

cf_pred <- c()
cc <- cf_matrix[8000001:nrow(cf_matrix),]
pred <- model %>% predict(cc)
cf_pred <- c(cf_pred,pred)
cf_pred <- as.data.frame(cf_pred)
cf_pred$len <- nchar(cf_seq$sequence)
write.table(cf_pred,'./data/cancer_data/bc_02_pred.txt',sep = '\t',quote = F,row.names = F)

####### get filter motif
 
seq_matrix <- cf_matrix[index,]
sequence <- cf_seq$sequence[index]

motif_value <- get_seqlogo(seq_matrix,sequence,model,0.7)
win_value <- motif_value[[2]]
motif <- motif_value[[1]]
pos_len <- motif_value[[3]]
len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))

add_NN <- function(aa)
{
  index <- which(nchar(aa)<19)
  for (i in index) {
    len <- nchar(aa[i])
    aa[i] <- paste0(aa[i],paste(rep('N',19-len),collapse = ''))
  }
  return(aa)
}
motif_nn <- lapply(1:length(motif), function(x) {motif[[x]][which(nchar(motif[[x]])==19)]})

filter_index <- which(len_m>10)
ggseqlogo(motif_nn[filter_index])

################ get depth

options (scipen = 999)

index <- which(cf_pred$cf_pred>=0.80)
all <- cbind(cf_seq[index,1:2],cf_pred[index,])
all$position <- gsub('>','',all$position)
all <- separate(all,'position',into = c('chrom','start','end'),sep='[:-]')
all$start <- as.integer(all$start)
all$end <- as.integer(all$end)
#all$mid_site <- round((all$start+all$end)/2)

mid_bed <- all[,c(1,2)]
colnames(mid_bed)[2] <- 'start'
mid_bed$end <- mid_bed$start+1000
mid_bed$start <- mid_bed$start-1000
mid_bed$point <- paste(mid_bed$chrom,paste(mid_bed$start,mid_bed$end,sep = '-'),sep = ':')
write.table(mid_bed[,4],'./data/test_2k.bed',sep = '\t',row.names = F,col.names = F,quote = F)
#system('bash /mnt/data1/qiting/cfDNA/TF_binding/data/get_depth.sh /mnt/data1/qiting/cfDNA/TF_binding/data/test_2k.bed /mnt/data1/qiting/cfDNA/TF_binding/data/cancer_data/picard/SRR14140136.sorted.markdup.bam')

mid_2k <- fread('./data/test_2k_depth.txt',sep = '\t')
tt <- matrix(mid_2k$V3,nrow = 2001,byrow = F)
#mid_2k <- apply(mid_bed,1,function(x) {cmd <- paste("samtools depth -a -r",paste(x[1],":",x[2],"-",x[3],sep = ''),"/mnt/data1/qiting/cfDNA/refer_expression/picard/SRR17478151.sorted.markdup.bam | awk '{print $3}'");system(cmd,intern = T)})
tt <- data.frame(mean_depth=rowMeans(tt))
tt$mean_depth <- tt$mean_depth/mean(tt$mean_depth)
tt$point <- seq(1,nrow(tt))
ggplot(tt, aes(x=point, y=mean_depth)) + geom_line(color="gray50", size=0.75)


######################## get position
refseq <- read.table('./data/gene_tss.txt',sep = '\t',header = T)
tss <- unique(refseq[,c(3,4,5,7,8,13)])

gene <- read.table('/mnt/data1/qiting/cfDNA/refer_expression/bed/gene_position_hg38.txt')
gene <- unique(gene[which(gene$V7=='protein_coding'),c(1,3,4,5,7,10)])

gene_tss <- merge(tss,gene,by.x='name2',by.y='V10')
gene_tss <- gene_tss[which(gene_tss$chrom %in% chromosomes),c(1,2,3,4,5,6)]
gene_tss$txStart <- gene_tss$txStart+1
#colnames(gene_tss) <- c('gene_id','chrom','strand','txStart','start','end')
gene_tss$tss <- gene_tss$txStart-2000
gene_tss$tss_end <- gene_tss$txStart+2000

temp <- gene_tss[-which(gene_tss$txStart==gene_tss$start),]
temp <- gene_tss


all[,7:9] <- NA
names(all)[c(5,7,8,9)] <- c('pred_value','up_TSS','down_TSS','ORF')
for (i in 1:nrow(all)) {
  chr <- all[i,1]
  start <- all[i,2]
  end <- all[i,3]
  a <- temp[which(temp$chrom==chr&temp$tss<start&temp$txStart>end),1]
  if(length(a)!=0){
    a <- unique(a)
    a <- paste(a,collapse  = '|')
    all[i,7] <- a
  }
  b <- temp[which(temp$chrom==chr&temp$txStart<start&temp$tss_end>end),1]
  if(length(b)!=0){
    b <- unique(b)
    b <- paste(b,collapse  = '|')
    all[i,8] <- b
  }
  c <- temp[which(temp$chrom==chr&temp$cdsStart<start&temp$cdsEnd>end),1]
  if(length(c)!=0){
    c <- unique(c)
    c <- paste(c,collapse  = '|')
    all[i,9] <- c
  }
  
}

View(table(all$up_TSS))
View(table(all$down_TSS))
View(table(all$ORF))

##############################################################
library(ChIPseeker)
library(ggplot2)
library(ggupset)
library(ggimage)
library(Guitar)
library(clusterProfiler)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb = TxDb.Hsapiens.UCSC.hg38.knownGene

index <- which(cf_pred$cf_pred>=0.995)
all_bed <- cf_bed[index,]
names(all_bed) <- c('chrom','start','end')
write.table(all_bed,'./result//bc_03_peak.bed',sep = '\t',quote = F,row.names = F,col.names = F)

stBedFiles <- list("/mnt/data1/qiting/cfDNA/TF_binding/result//bc_03_peak.bed")
p <- GuitarPlot(txTxdb = txdb, 
                stBedFiles = stBedFiles, 
                headOrtail = FALSE,
                enableCI = FALSE, 
                mapFilterTranscript = TRUE, 
                pltTxType = c("mrna"), 
                stGroupName = "KO1")

p <- p + ggtitle(label = "Distribution on mRNA") + 
  theme(plot.title = element_text(hjust = 0.5)) + theme_bw()

f <- GRanges(seqnames = all_bed$chrom,IRanges(start = all_bed$start,end = all_bed$end))
covplot(f)

x = annotatePeak(f, tssRegion=c(-3000, 3000), TxDb=txdb,ignoreUpstream = T,ignoreDownstream = T)
plotAnnoPie(x)
plotAnnoBar(x)
plotAnnoPie.csAnno(x)
plotDistToTSS(x)
upsetplot(x, vennpie=TRUE)


peakHeatmap(f,TxDb=txdb, upstream=3000, downstream=3000)


require(clusterProfiler)
bedfile=getSampleFiles()
# 将bed文件读入（readPeakFile是利用read.delim读取，然后转为GRanges对象）
seq=lapply(stBedFiles, readPeakFile)

genes=lapply(seq, function(i) seq2gene(i, c(-1000, 3000), 3000, TxDb=txdb))
go <- enrichGO(gene=genes[[1]],OrgDb='org.Hs.eg.db',ont= "CC",pAdjustMethod="BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
dotplot(go)

cc = compareCluster(geneClusters = genes, fun="enrichKEGG", organism="hsa")
dotplot(cc, showCategory=10)

