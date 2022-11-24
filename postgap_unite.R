
#rs <- list.files(path = "/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/rsids/")
#ch <- c(15:22)
#for (k in ch) {
#  for (i in 1:length(rs)) {
#    if (i == 1) {
#      init <- read.table(paste('/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/rsids/',rs[i],sep=''),sep= '\t',header = T)
#    }
#    else {
#      print(c(rs[i],k,i))
#      f <- try(read.table(paste('/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/rsids/',rs[i],sep=''),sep= '\t',header = T))
#      if (inherits(f, 'try-error')) next
#      if (nrow(f) != 0 ) {
#            if (f$chrom[1] == k) {
#                init <- rbind(init,f)
#            }
#        }
#      }
#    
#  }
#  write.table(init,paste("/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_",k,"_postgap.tsv",sep=''), sep='\t',quote = F,row.names = F,col.names = T)
#}

ch <- c(1:22)
for (i in 1:length(ch)) {
    if (i == 1) {
      init <- read.table(paste('/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_',ch[i],'_postgap.tsv',sep=''),sep= '\t',header = T)
    }
    else {
      print(c(i))
      f <- try(read.table(paste('/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_',ch[i],'_postgap.tsv',sep=''),sep= '\t',header = T))
      if (inherits(f, 'try-error')) next
      init <- rbind(init,f)
        }
}
init <- unique(init)
write.table(init,paste("/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_postgap.tsv",sep=''), sep='\t',quote = F,row.names = F,col.names = T)


EA <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_postgap.tsv",sep = '\t', header = T)
colnames(EA)[1] <- 'ld_snp_rsID'
EA_test <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/20116_0.gwas.imputed_v3.both_sexes_FILTER.tsv",sep = '\t', header = T)

EA_test <- EA_test[,c(3,12,9)]
colnames(EA_test) <- c('ld_snp_rsID','p.value','beta')
EA <- as.data.frame(merge(EA,EA_test,by='ld_snp_rsID'))
EA_test <- EA[,c(1,83,84)]
EA_1 <- EA[,c(1:25)]
EA_1$disease_name <- 'Manual'
EA_1$disease_efo_id <- 'EFO_Manual'
EA_1$score <- EA$score
EA_1$rank <- 1
EA_1$r2 <- 1
EA_1$cluster_id <- 6.715048e+18
EA_1$gwas_source <- 'Manual'
EA_1$gwas_snp <- EA_test$ld_snp_rsID
EA_1$gwas_pvalue <- EA_test$p.value
EA_1$gwas_pvalue_description <- 'Manual'
EA_1$gwas_odds_ratio <- 'None'
EA_1$gwas_odds_ratio_ci_start <- 'None'
EA_1$gwas_odds_ratio_ci_end <- 'None'
EA_1$gwas_beta <- EA_test$beta
EA_1$gwas_size <- 1000
EA_1$gwas_pmid <- 'PMID000'
EA_1$gwas_study <- 'Manual'
EA_1$gwas_reported_trait <- 'Manual'
EA_1$ld_snp_is_gwas_snp <- 1

EA_2 <- EA[,c(27:82)]
EA_1 <- cbind(EA_1,EA_2)
colnames(EA_1)[which(names(EA_1) == "GTEx_Cells_EBV.transformed_lymphocytes")] <- "GTEx_Cells_EBV-transformed_lymphocytes"


write.table(EA_1,paste("/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_postgap.tsv",sep=''), sep='\t',quote = F,row.names = F,col.names = T)

