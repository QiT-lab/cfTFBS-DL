##############################################
#### calculate AUC and plot ROC curve
##############################################
library(pROC)
library(ggplot2)

oneHotEncode <- function(bio_seq)
{
  test <- data.frame(row.names =c('A','T','C','G'), col1=c(1,0,0,0),col2=c(0,1,0,0),col3=c(0,0,1,0),col4=c(0,0,0,1))
  bio_char <- strsplit(bio_seq,'')[[1]]
  df <- unlist(lapply(1:length(bio_char), function(i) {return(test[bio_char[i],])}))
  return(df)
}


cnn_model <- function(shuffle_seq,k)
{
  pos_list <- lapply(shuffle_seq$sequence, oneHotEncode)
  pos_matrix <- do.call(rbind,lapply(1:length(pos_list), function(x) {get_samelength(pos_list[[x]])}))
  
  neg_list <- lapply(shuffle_seq$neg_sequence, oneHotEncode)
  neg_matrix <- do.call(rbind,lapply(1:length(neg_list), function(x) {get_samelength(neg_list[[x]])}))
  
  set.seed(12)
  q <- round(nrow(shuffle_seq)/10)
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
    layer_conv_1d(filters =128, kernel_size = k, activation = "relu", input_shape = c(276,1),strides = 4) %>%
    layer_max_pooling_1d(pool_size = 2) %>%
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
  
  history <- model %>% fit(train_data,train_label,epochs = 20,batch_size = 32,validation_split = 0.1)
  result <- model %>% evaluate(test_data,test_label)
  print(result)
  pred <- model %>% predict(test_data)
  res <- roc(test_label,as.vector(pred),aur=TRUE,ci=TRUE)
  #save_model_hdf5(model,'./model/CTCF_breast_model.hdf5')
  return(list(model,res))
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


gen_neg_seq <- function(chromosomes,chrom_size,temp,width,seq)
{
  options(scipen = 99999)
  gc <- temp$gc
  random_gc=0
  neg_seq <- 'NNNN'
  width <- width - 1 
  while(length(grep('N',neg_seq))!=0 | abs(random_gc-gc) > 0.1 ){
    random_chrom <- sample(chromosomes,1)
    chrom_length <- chrom_size[which(chrom_size$V1==random_chrom),2]
    random_start <- sample(1:(chrom_length-width),1)
    random_end <- random_start+width
    same_chr <- seq[which(seq$chrom==random_chrom),]
    t <- apply(same_chr, 1, function(x) {intersect(random_start:random_end,x[2]:x[3])})
    if(length(t)==0){ 
      cmd <- paste('samtools faidx /mnt/public/reference/human/hg38/bwa/hg38.fa',paste0(random_chrom,':',random_start,'-',random_end),sep = ' ')
      neg_seq <- system(cmd,intern = TRUE)
      position <- neg_seq[1]
      if(width>60){
        neg_seq <- toupper(paste0(neg_seq[2],neg_seq[3]))
      }else{
        neg_seq <- toupper(neg_seq[2])
      }
      random_gc <- calculate_gc(neg_seq)
    }
  }
  temp$neg_position <- position
  temp$neg_sequence <- neg_seq
  temp$neg_gc <- random_gc
  return(temp)
}

get_shuffle_seq <- function(bed,gen_neg_seq)
{
  seq <- getfasta(bed)
  seq <- seq[which(!grepl('N',seq$sequence)),]
  
  seq$gc <- unlist(lapply(1:nrow(seq),function(i) {calculate_gc(seq$sequence[i])}))
  seq$position <- gsub('>','',seq$position)
  seq <- separate(seq,'position',into = c('chrom','start','end'),sep='[:-]')
  seq$start <- as.integer(seq$start)
  seq$end <- as.integer(seq$end)
  seq$width <- seq$end-seq$start
  
  chromosomes <- c(paste0("chr", 1:22),'chrX','chrY')
  chrom_size <- read.table('./data/chromsize.txt')
  chrom_size <- chrom_size[which(chrom_size$V1 %in% chromosomes),]
  
  
  library(parallel)
  detectCores()  
  clus <- makeCluster(80)
  clusterExport(clus,deparse(substitute(gen_neg_seq)))
  clusterExport(clus,"chromosomes",envir = environment())
  clusterExport(clus,"chrom_size",envir = environment())
  clusterExport(clus,"seq",envir = environment())
  #clusterExport(clus,"width",envir = environment())
  clusterExport(clus,deparse(substitute(calculate_gc)))
  
  shuffle_seq <- parLapply(clus,1:nrow(seq), function(i) gen_neg_seq(chromosomes,chrom_size,seq[i,],seq[i,6],seq))
  shuffle_seq <- do.call(rbind,shuffle_seq)
  stopCluster(clus)
  return(shuffle_seq)
}

shuffle_seq <- import_fst('./data/gtrd_cell_line/CTCF_A549_shuffle_seq.fst')
all <- cnn_model(shuffle_seq,k=76)
model1 <- all[[1]]
save_model_hdf5(model1,'./model/CTCF_A549_model.hdf5')
res1 <- all[[2]]
print(res1)


shuffle_seq <- import_fst('./data/gtrd_cell_line/CTCF_MCF7_shuffle_seq.fst')
all <- cnn_model(shuffle_seq,k=76)
model2 <- all[[1]]
save_model_hdf5(model2,'./model/CTCF_MCF7_model.hdf5')
res2 <- all[[2]]
print(res2)

shuffle_seq <- import_fst('./data/gtrd_cell_line/CTCF_monocytes_shuffle_seq.fst')
all <- cnn_model(shuffle_seq,k=76)
model3 <- all[[1]]
save_model_hdf5(model3,'./model/FOXA1_breast_model.hdf5')
res3 <- all[[2]]


shuffle_seq <- import_fst('./data/gtrd_cell_line/FOXA1_MCF7_shuffle_seq.fst')
all <- cnn_model(shuffle_seq,k=44)
model4 <- all[[1]]
save_model_hdf5(model4,'./model/FOXA1_MCF7_model.hdf5')
res4 <- all[[2]]


shuffle_seq <- import_fst('./data/gtrd_cell_line/FOXA1_A549_shuffle_seq.fst')
all <- cnn_model(shuffle_seq,k=44)
model5 <- all[[1]]
save_model_hdf5(model5,'./model/FOXA1_A549_model.hdf5')
res5 <- all[[2]]


#res <- list(`CTCF lung: AUC=0.97`=res1,`CTCF breast: AUC=0.91`=res2, `FOXA1 breast: AUC=0.85`=res3,`CTCF monocyte: AUC=0.90`=res4,`CTCF lymphocyte: AUC=0.88`=res5)
res <- list(`CTCF lung: AUC=0.82`=res1,`CTCF breast: AUC=0.84`=res2, `FOXA1 breast: AUC=0.83`=res4,`FOXA1 lung: AUC=0.84`=res5)
ggroc(res,legacy.axes = T)+
  #geom_abline(intercept = 0,slope=270)+
  annotate(geom = 'segment',x=0,y=0,xend = 1,yend = 1,lty=2)+
  #annotate(geom = 'text',label='AUC=0.8584',size=5,color='red',x=0.3,y=0.5)+
  my_theme+
  theme(panel.border = element_rect(size=1,fill = 'transparent'),
        legend.position = c(0.8,0.3))

ggsave('./result/five_sample_roc.pdf',height = 4, width = 6)
###########################################
options(scipen = 99999)
setwd('/mnt/data1/qiting/cfDNA/TF_binding/')
dir <- '/mnt/data1/qiting/cfDNA/TF_binding/data/'
peak <- fread('./data/gtrd_cell_line/human_meta_clusters.interval.gz')
names(peak)[1] <- 'CHROM'
chromosomes <- c(paste0("chr", 1:22),'chrX','chrY')
peak <- peak %>% filter_dt(CHROM %in% chromosomes)
tf <- c('CTCF [11]',"FOXA1 (HNF-3&#945;)")
cell_list <- list()
cell_list[[2]] <- c('MCF7 \\(Invasive ductal breast carcinoma\\)','A549 \\(lung carcinoma\\)')
cell_list[[1]] <- c('MCF7 \\(Invasive ductal breast carcinoma\\)','A549 \\(lung carcinoma\\)','blood lymphocytes','blood monocytes')
for (f in 1:2) {
  print(tf[f])
  data <- peak %>% filter_dt(tfTitle==tf[f])
  #cellset <- unique(unlist(strsplit(data$cell.set,';')))
  cell <- cell_list[[f]]
  for (i in 3:length(cell)) {
    print(cell[i])
    cell_line <- cell[i]
    tt <- data[grep(cell_line,data$cell.set),]
    tt$len <- tt$END-tt$START
    temp <- tt[which(tt$len<70&tt$len>40),]
    bed <- unique(temp[,c(1,2,3)])
    shuffle_seq <- get_shuffle_seq(bed,gen_neg_seq)
    tf_id <- strsplit(tf[f],' ')[[1]][1]
    cell_id <- strsplit(cell_line,' ')[[1]][1]
    if(cell_id=='blood'){
      export_fst(shuffle_seq,paste0('./data/gtrd_cell_line/',tf_id,gsub('blood ','_',cell_line),'_shuffle_seq.fst'))
    }else{
      export_fst(shuffle_seq,paste0('./data/gtrd_cell_line/',tf_id,'_',cell_id,'_shuffle_seq.fst'))
    }
  }
}


