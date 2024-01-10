######################################################
######### Visualizing Convolutional Layer
######### get motif from filter
######################################################

library(ggseqlogo)


get_seqlogo <- function(seq_matrix, sequence,model,p,n)
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
    tt[[k]] <- apply(aa,1,function(x) {substr(sequence[x[1]],x[2],x[2]+n)})
  }
  return(list(tt,win_value,pos_len))
}


shuffle_seq <- import_fst('./data/gtrd_cell_line/CTCF_MCF7_shuffle_seq.fst')
model1 <- load_model_hdf5('./model/CTCF_MCF7_model.hdf5')
pos_list <- lapply(shuffle_seq$sequence, oneHotEncode)
pos_matrix <- do.call(rbind,lapply(1:length(pos_list), function(x) {get_samelength(pos_list[[x]])}))
motif_value <- get_seqlogo(pos_matrix,shuffle_seq$sequence,model1,p=1,n=18)
motif <- motif_value[[1]]
cf_motif <- lapply(1:length(motif), function(x) {motif[[x]][which(nchar(motif[[x]])==19)]})
len_m <- unlist(lapply(1:length(cf_motif), function(x) {length(cf_motif[[x]])}))
#len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))
filter_index <- which(len_m>1000)
ggseqlogo(cf_motif[filter_index])
ggseqlogo(cf_motif[64])
ggsave('./result/CTCF_MCF7_motif.pdf',width = 6,height = 4)


shuffle_seq <- import_fst('./data/gtrd_cell_line/CTCF_A549_shuffle_seq.fst')
model2 <- load_model_hdf5('./model/CTCF_A549_model.hdf5')
pos_list <- lapply(shuffle_seq$sequence, oneHotEncode)
pos_matrix <- do.call(rbind,lapply(1:length(pos_list), function(x) {get_samelength(pos_list[[x]])}))
motif_value <- get_seqlogo(pos_matrix,shuffle_seq$sequence,model2,p=1,n=18)
motif <- motif_value[[1]]
cf_motif <- lapply(1:length(motif), function(x) {motif[[x]][which(nchar(motif[[x]])==19)]})
len_m <- unlist(lapply(1:length(cf_motif), function(x) {length(cf_motif[[x]])}))
#len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))
filter_index <- which(len_m>100)
ggseqlogo(cf_motif[filter_index])
ggseqlogo(cf_motif[filter_index[13]])
ggsave('./result/CTCF_A549_motif.pdf',width = 6,height = 4)


shuffle_seq <- import_fst('./data/gtrd_cell_line/FOXA1_MCF7_shuffle_seq.fst')
model3 <- load_model_hdf5('./model/FOXA1_MCF7_model.hdf5')
pos_list <- lapply(shuffle_seq$sequence, oneHotEncode)
pos_matrix <- do.call(rbind,lapply(1:length(pos_list), function(x) {get_samelength(pos_list[[x]])}))
motif_value <- get_seqlogo(pos_matrix,shuffle_seq$sequence,model3,p=1,n=10)
motif <- motif_value[[1]]
cf_motif <- lapply(1:length(motif), function(x) {motif[[x]][which(nchar(motif[[x]])==11)]})
len_m <- unlist(lapply(1:length(cf_motif), function(x) {length(cf_motif[[x]])}))
#len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))
filter_index <- which(len_m>1000)
ggseqlogo(cf_motif[filter_index])
ggseqlogo(cf_motif[filter_index[3]])
ggsave('./result/FOXA1_MCF7_motif.pdf',width = 6,height = 4)
ggseqlogo(cf_motif[filter_index[36]])
ggsave('./result/FOXA1_MCF7_motif_2.pdf',width = 6,height = 4)

shuffle_seq <- import_fst('./data/gtrd_cell_line/FOXA1_A549_shuffle_seq.fst')
model4 <- load_model_hdf5('./model/FOXA1_A549_model.hdf5')
pos_list <- lapply(shuffle_seq$sequence, oneHotEncode)
pos_matrix <- do.call(rbind,lapply(1:length(pos_list), function(x) {get_samelength(pos_list[[x]])}))
motif_value <- get_seqlogo(pos_matrix,shuffle_seq$sequence,model4,p=0.7,n=10)
motif <- motif_value[[1]]
cf_motif <- lapply(1:length(motif), function(x) {motif[[x]][which(nchar(motif[[x]])==11)]})
len_m <- unlist(lapply(1:length(cf_motif), function(x) {length(cf_motif[[x]])}))
#len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))
filter_index <- which(len_m>20)
ggseqlogo(cf_motif[filter_index])
ggseqlogo(cf_motif[filter_index[35]])
ggsave('./result/FOXA1_A549_motif.pdf',width = 6,height = 4)


shuffle_seq <- read.table('./data/CTCF_data/GSE46832.CTCF.lymphocyte_blood_shuffle_81.txt',sep = '',header = T)
model5 <- load_model_hdf5('./model/CTCF_lymphocyte_model.hdf5')
pos_matrix <- do.call(rbind,lapply(shuffle_seq$sequence, oneHotEncode))
motif_value <- get_seqlogo(pos_matrix,shuffle_seq$sequence,model5,p=0.9,n=18)
motif <- motif_value[[1]]
len_m <- unlist(lapply(1:length(motif), function(x) {length(motif[[x]])}))
filter_index <- which(len_m>100)
ggseqlogo(motif[filter_index])
ggseqlogo(motif[filter_index[7]])
ggsave('./result/CTCF_lymphocyte_motif.pdf',width = 6,height = 4)
