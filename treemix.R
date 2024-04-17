library(dplyr)
library(tidyr)

AC <- 'path with ACs for each population'
input <- 'input treemix file'

for (i in c('cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')) {
  if (i == 'cl1') {

    test <- read.table(paste(AC,i,"_all.txt",sep=''),sep = '\t', header = T)
    test <- test %>% unite('locus',c('locus','alleles'),sep='_')
    colnames(test)[2] <- 'AC'
    test$AC <- gsub('\\[|\\]','',test$AC)
    
  }
  else {
    test_ <- read.table(paste(AC,i,"_all.txt",sep=''),sep = '\t', header = T)
    test_ <- test_ %>% unite('locus',c('locus','alleles'),sep='_')
    colnames(test_)[2] <- 'AC'
    test_$AC <- gsub('\\[|\\]','',test_$AC)
    test <- as.data.frame(merge(test,test_,by='locus'))
  }
}

colnames(test) <- c('locus','cl1','cl2','cl3','cl4','cl5','cl6','ASW','ACB','LWK','GWD','YRI','MSL','ESN','PUR','CLM','MXL','PEL','PJL','GIH','BEB','ITU','STU','CHB','JPT','CHS','CDX','KHV','CEU','GBR','FIN','IBS','TSI')

write.table(test[,2:ncol(test)],input,sep=' ',quote = F,row.names = F,col.names = T)

