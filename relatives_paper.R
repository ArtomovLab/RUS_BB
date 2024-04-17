library(reshape2)
library(dplyr)
library(tidyr)
library(igraph)

duplicates <- 'txt file with IDs of duplicated samples'
king_matrix <- 'plink kinship matrix'
king_matrix_id <- 'IDs of plink kinship matrix'
esse_pheno <- 'txt file of Russian phenotypes'

duplicates <-  read.table(duplicates,sep = '\t', header = F)

king <-  read.table(king_matrix,sep = '\t', header = F)
king_id <-  read.table(king_matrix_id,sep = '\t', header = F)
colnames(king) <- king_id$V1
row.names(king) <- king_id$V1
king$sample = king_id$V1

tmp_king <- melt(king,id=ncol(king)) 
tmp_king <- tmp_king %>% filter(value != 0)
tmp_king <- tmp_king %>% filter(sample != variable) %>% dplyr::arrange(desc(value))
tmp_king <- tmp_king %>% filter(sample %!in% duplicates$V1 & variable %!in% duplicates$V1)
#write.table(tmp_king,"plink2.king_long", sep='\t',quote = F,row.names = F,col.names = F)

tmp_king_relatives <- tmp_king %>% filter(value >= 0.08838835)
test <- tmp_king %>% filter(value <= 0.08838835 & value >= 0.04419417)
tmp_king_relatives <- tmp_king_relatives %>% mutate(gen = case_when(value > 0.04419417 & value < 0.08838835 ~ '3rd',
                                                                    value > 0.08838835 & value < 0.1767767 ~ '2nd',
                                                                    value > 0.1767767 & value < 0.3535534 ~ '1st',
                                                                    value > 0.3535534 ~ 'twin'))
table(tmp_king_relatives$gen)

colnames(tmp_king_relatives)[3] <- 'weight'
tmp_king_relatives$weight <- round(tmp_king_relatives$weight,3)
colnames(tmp_king_relatives)[1:2] <- c('from','to')
##
General <- read.table(esse_pheno,sep = '\t', header = T)

ddd <- General %>% 
  dplyr::select(NP,BMI,AGE,SBP,DBP,HEI,WEI,DATEBIRTH,DateInform,PATIENTFIO,SEX)
dd <- Gen_gwas %>% dplyr::select(ID,AGE,SEX)
colnames(dd)[1] <- 'from'
tmp_king_relatives <- as.data.frame(merge(tmp_king_relatives,dd,by='from',all.x=T))
tmp_king_relatives <- tmp_king_relatives %>% unite('from',c('from','AGE','SEX'),sep='_')

colnames(dd)[1] <- 'to'
tmp_king_relatives <- as.data.frame(merge(tmp_king_relatives,dd,by='to',all.x=T))
tmp_king_relatives <- tmp_king_relatives %>% unite('to',c('to','AGE','SEX'),sep='_')


g <- graph_from_data_frame(tmp_king_relatives, directed=FALSE)
#g_twin <- graph_from_data_frame(tmp_king_relatives, directed=FALSE)
#plot(g_twin,vertex.size=0.1)

E(g)$type <- tmp_king_relatives$gen
E(g)$width <- 2
dg <- decompose.graph(g) 

length(dg)

L <- 0
I <- 0
for (i in 1:length(dg)) {
  l <-  length(E(dg[[i]]))
  print(l)
  if (l > L) {
    L <- l
    I <- i
  }
}
print(I)


n <- length(dg)
for (j in 1:(n-1)) {
  for(i in 1:(n-j)) {
    if (length(E(dg[[i]])) < length(E(dg[[i+1]]))  ) {
      temp <- dg[[i+1]]
      dg[[i+1]] <- dg[[i]]
      dg[[i]] <- temp
    }
  }
}


for (i in 1:length(dg)) {
  col <- gsub('twin',1,E(dg[[i]])$type)
  col <- gsub('1st',2,col)
  col <- gsub('2nd',3,col)
  col <- gsub('3rd',4,col)
  col <- as.numeric(col)
  png(paste(i,'.png',sep=''))
  plot(dg[[i]],edge.color=c("red", "blue",'green','yellow')[col]) 
  dev.off()
}

