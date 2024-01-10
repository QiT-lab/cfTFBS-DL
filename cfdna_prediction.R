#######################################################
################## predict TFBS in cfDNA sequence
################ get TFBS motif
################ get depth in TFBS motif
########################################################
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

oneHotEncode <- function(bio_seq)
{
  test <- data.frame(row.names =c('A','T','C','G'), col1=c(1,0,0,0),col2=c(0,1,0,0),col3=c(0,0,1,0),col4=c(0,0,0,1))
  bio_char <- strsplit(bio_seq,'')[[1]]
  df <- unlist(lapply(1:length(bio_char), function(i) {return(test[bio_char[i],])}))
  return(df)
}

get_samelength <- function(aa)
{
  ll <- length(aa)
  if(ll < 276){
    aa <- c(aa,rep(0,(276-ll)))
  }
  return(aa)
}

get_cf_matrix <- function(temp_seq)
{
  library(parallel)
  detectCores()  
  clus <- makeCluster(80)
  clusterExport(clus,deparse(substitute(oneHotEncode)))
  #clusterExport(clus,"cf_seq",envir = environment())
  cf_list <- parLapply(clus,temp_seq$sequence, oneHotEncode)
  stopCluster(clus)
  
  cf_matrix <- do.call(rbind,lapply(1:length(cf_list), function(x) {get_samelength(cf_list[[x]])}))
  #cf_matrix[which(is.na(cf_matrix),arr.ind = T)] <- 0
  return(cf_matrix)
}

library(tidyfst)
options (scipen = 9999)
setwd('/mnt/data1/qiting/cfDNA/TF_binding/')
dir <- '/mnt/data1/qiting/cfDNA/TF_binding/data/'
chromosomes <- c(paste0("chr", 1:22),'chrX','chrY')
file <- c('bc_02.bed','bc_03.bed','lung_G2.bed','lung_A2.bed','healthy_F02.bed')
# for (f in file[2:4]) {
#   print(f)
#   cfdna <- data.table::fread(paste0('./data/cancer_data/bed/',f))
#   cf_bed <- unique(cfdna[which(cfdna$V4<82),c(1,2,3)])
#   cf_bed <- cf_bed[which((cf_bed$V1 %in% chromosomes)&(cf_bed$V2>0)),]
#   cf_seq <- getfasta(cf_bed)
#   cf_seq <- cf_seq[-which(grepl('N',cf_seq$sequence)),]
#   write.table(cf_seq,paste0('./result/seq/',sub("*.bed", "", f),'_cf_seq.txt'),quote = F,sep = '\t',row.names = F)
#   
#   aa <- c(seq(1,nrow(cf_seq),1000000),nrow(cf_seq)+1)
#   cf_matrix <- matrix(nrow = 0, ncol = 324)
#   for (i in 1:(length(aa)-1)) {
#     print(i)
#     temp_seq <- cf_seq[aa[i]:(aa[i+1]-1),]
#     temp_matrix <- get_cf_matrix(temp_seq)
#     cf_matrix <- rbind(cf_matrix,temp_matrix)
#   }
#   write.table(cf_matrix,paste0('./result/matrix/',sub("*.bed", "", f),'_81_matrix.txt'),quote = F,sep = '\t',row.names = F,col.names = F)
# }


for (f in file[4]) {
  print(f)
  cfdna <- data.table::fread(paste0('./data/cancer_data/bed/',f))
  cf_bed <- cfdna %>% filter_dt(V4 < 70 & V4 > 40) %>% select_dt(c(1,2,3)) %>% unique()
  cf_bed <- cf_bed %>% filter_dt((V1 %in% chromosomes) & (V2 > 0))
  cf_seq <- getfasta(cf_bed)
  cf_seq <- cf_seq %>% filter_dt(!grepl('N',cf_seq$sequence))
  export_fst(cf_seq,paste0('./result/seq/',sub("*.bed", "", f),'_cf_seq_69.fst'))
  #cf_seq <- import_fst('./result/seq/bc_02_cf_seq.fst')
  
  aa <- c(seq(1,nrow(cf_seq),1000000),nrow(cf_seq)+1)
  cf_matrix <- matrix(nrow = 0, ncol = 276)
  for (i in 1:(length(aa)-1)) {
    print(i)
    temp_seq <- cf_seq[aa[i]:(aa[i+1]-1),]
    temp_matrix <- get_cf_matrix(temp_seq)
    cf_matrix <- rbind(cf_matrix,temp_matrix)
  }
  cf_matrix <- as.data.frame(cf_matrix)
  export_fst(cf_matrix,paste0('./result/matrix/',sub("*.bed", "", f),'_69_matrix.fst'))
}




model <- load_model_hdf5('./model/FOXA1_A549_model.hdf5')
cf_matrix <- import_fst('./result/matrix/lung_G2_69_matrix.fst')
cf_matrix <- as.matrix(cf_matrix)
cf_seq <- import_fst('./result/seq/lung_G2_cf_seq_69.fst')
 
bb <- c(seq(1,nrow(cf_matrix),2000000),nrow(cf_matrix)+1)
pred <- c()
for (j in 1:(length(bb)-1)) {
  temp_pred <- model %>% predict(cf_matrix[bb[j]:(bb[j+1]-1),])
  pred <- c(pred,temp_pred)
}
cf_pred <- cbind(as.data.frame(pred),nchar(cf_seq$sequence))
names(cf_pred) <- c('score','len')
export_fst(cf_pred,'./result/prediction/FOXA1_A549_lung_G2_predict.fst')
#write.table(cf_pred,'./result/prediction/CTCF_MCF7_bc_02_predict.txt',sep = '\t',quote = F,row.names = F)

#################################### ##################
######the relation of prediction score and cfDNA length
#######################################################
cf_pred <- read.table('./result/prediction/CTCF_breast_bc_02_predict.txt',sep = '\t',header = T)
temp <- cf_pred[which(cf_pred$score>0.8),]
mean_score <- aggregate(score ~ len, data = temp, FUN = mean)
ggplot(mean_score, aes(x = len, y = score)) +   
  geom_point(colour = "red", size = 2) +  
  geom_smooth(method = "loess") +
  #ggtitle("My Line Chart") +  
  xlab("cfDNA fragment length") +  
  ylab("mean prediction score") + 
  my_theme
ggsave('./result/CTCF_lung_lung_G2_mean_score_len.pdf',width = 6,height = 4)


############################### get cfdna filter
library(ggseqlogo)
cf_pred <- read.table('./result//FOXA1_breast_bc_02_predict.txt',sep = '\t',header = T)
index <- which(cf_pred$score>=0.98)
seq_matrix <- cf_matrix[index,]
sequence <- cf_seq$sequence[index]

motif_value <- get_seqlogo(seq_matrix,sequence,model,p=0.9,n=10)
motif <- motif_value[[1]]
cf_motif <- lapply(1:length(motif), function(x) {motif[[x]][which(nchar(motif[[x]])==11)]})
len_m <- unlist(lapply(1:length(cf_motif), function(x) {length(cf_motif[[x]])}))
filter_index <- which(len_m>1000)
ggseqlogo(cf_motif[filter_index])
ggseqlogo(cf_motif[filter_index[7]])
ggsave('./result/CTCF_lung_motif.pdf',width = 6,height = 4)


###################################  get depth
options (scipen = 999)

cf_pred <- read.table('./result/prediction/CTCF_lung_lung_G2_predict.txt',sep = '\t',header = T)
cf_seq <- import_fst('./result/seq/lung_G2_cf_seq.fst')
index <- which(cf_pred$score>=0.98)
all <- cbind(cf_seq[index,1:2],cf_pred[index,])
all$position <- gsub('>','',all$position)
all <- separate(all,'position',into = c('chrom','start','end'),sep='[:-]')
all$start <- as.integer(all$start)
all$end <- as.integer(all$end)
all$start <- floor((all$start+all$end)/2)

mid_bed <- all[,c(1,2)]
colnames(mid_bed)[2] <- 'start'
mid_bed$end <- mid_bed$start+1000
mid_bed$start <- mid_bed$start-1000
mid_bed$point <- paste(mid_bed$chrom,paste(mid_bed$start,mid_bed$end,sep = '-'),sep = ':')
write.table(mid_bed[,4],'./data/test_2k.bed',sep = '\t',row.names = F,col.names = F,quote = F)

system('/bin/bash /mnt/data1/qiting/cfDNA/TF_binding/data/get_depth.sh /mnt/data1/qiting/cfDNA/TF_binding/data/test_2k.bed /mnt/data1/qiting/cfDNA/TF_binding/data/cancer_data/bam/bc_02_sorted.bam')

mid_2k <- fread('./data/test_2k_depth.txt',sep = '\t')
tt <- matrix(mid_2k$V3,nrow = 2001,byrow = F)
tt <- data.frame(mean_depth=rowMeans(tt))
tt$mean_depth <- tt$mean_depth/mean(tt$mean_depth)
tt$point <- seq(-1000,1000)
ggplot(tt, aes(x=point, y=mean_depth)) + 
  xlab("cfDNA upstream and downstream") +  
  ylab("mean depth") + 
  geom_line(color="black", size=0.5) +
  my_theme
ggsave('./result/depth/CTCF_lung_lung_G2_depth.pdf',width = 6,height = 4)


