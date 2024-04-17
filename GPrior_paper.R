library(tidyr)
library(dplyr)
library(ggplot2)

#SMOKING never
SMOKING <- 'list of genes from GWAS catalog for Smoking Initiation phenotype'
SMNE_train <- 'list of genes to train GPrior Smoking Initiation'
SMNE_gprior <- 'GPrior results for Smoking Initiation phenotype'
SMOKING_ever <- 'list of genes from GWAS catalog for Smoking Cessation phenotype'
SmokingCessation_train <- 'list of genes to train GPrior Smoking Cessation'
SmokingCessation_train_gprior <- 'GPrior results for Smoking Cessation phenotype'

SMOKING <- read.table(paste(SMOKING,sep=''),sep = '\t', header = T,quote = '"')
SMOKING <- SMOKING %>% separate('P.value', c('P.value1','P.value2'), sep= 'x')
SMOKING <- SMOKING %>% separate('P.value2', c('P.value2','P.value3'), sep= '-')
SMOKING$P.value1 <- as.numeric(SMOKING$P.value1)
SMOKING$P.value2 <- as.numeric(SMOKING$P.value2)
SMOKING$P.value3 <- as.numeric(SMOKING$P.value3)
SMOKING <- SMOKING %>% mutate(P.value = P.value1*(P.value2^(-P.value3)))
SMOKING <- SMOKING[!grepl(',',SMOKING$Mapped.gene),]
SMOKING <- SMOKING %>% filter(Mapped.gene != "'-")
SMOKING <- SMOKING %>% filter(P.value < 1e-6)
#length(unique(SMOKING$Mapped.gene))
#SMOKING <- SMOKING %>% filter(Mapped.gene %in% GPrior$gene_symbol)
SMOKING <- SMOKING %>% filter(Reported.trait %in% c('Smoking initiation (ever regular vs never regular)'))
SMOKING <- as.data.frame(unique(SMOKING$Mapped.gene))
colnames(SMOKING) <- 'gene_symbol'
write.table(SMOKING,paste(SMNE_train),sep='\t',quote = F,row.names = F)

GPrior_res <- read.table(paste(SMNE_gprior,sep=''),sep = '\t', header = T)
df2 <- GPrior_res %>%  gather(key,value,2:7) %>%
  filter(gene_symbol == 'PTK2')

GPrior_res %>% gather(key,value,2:7) %>%
  filter(value != -1) %>% 
  ggplot(aes(value))+
  geom_histogram()+
  facet_wrap(~key)+
  theme_classic()+
  geom_vline(data = df2, mapping = aes(xintercept = value)) 

##SMOKING Ever
SMOKING <- read.table(paste(SMOKING_ever,sep=''),sep = ',', header = T,quote = '"')
SMOKING <- SMOKING %>% separate('P.value', c('P.value1','P.value2'), sep= 'x')
SMOKING <- SMOKING %>% separate('P.value2', c('P.value2','P.value3'), sep= '-')
SMOKING$P.value1 <- as.numeric(SMOKING$P.value1)
SMOKING$P.value2 <- as.numeric(SMOKING$P.value2)
SMOKING$P.value3 <- as.numeric(SMOKING$P.value3)
SMOKING <- SMOKING %>% mutate(P.value = P.value1*(P.value2^(-P.value3)))
SMOKING <- SMOKING[!grepl(',',SMOKING$Mapped.gene),]
SMOKING <- SMOKING %>% filter(Mapped.gene != "'-")
SMOKING <- SMOKING %>% filter(P.value < 1e-6)
#length(unique(SMOKING$Mapped.gene))
#SMOKING <- SMOKING %>% filter(Mapped.gene %in% GPrior$gene_symbol)

SMOKING <- SMOKING %>% filter(Reported.trait %in% c('Smoking cessation'))
SMOKING <- as.data.frame(unique(SMOKING$Mapped.gene))
colnames(SMOKING) <- 'gene_symbol'
write.table(SMOKING,paste(SmokingCessation_train),sep='\t',quote = F,row.names = F)

GPrior_res <- read.table(paste(SmokingCessation_train_gprior,sep=''),sep = '\t', header = T)
df2 <- GPrior_res %>%  gather(key,value,2:7) %>%
  filter(gene_symbol == 'ACSM4')

GPrior_res %>% gather(key,value,2:7) %>%
  filter(value != -1) %>%
  ggplot(aes(value))+
  geom_histogram()+
  facet_wrap(~key)+
  theme_classic()+
  geom_vline(data = df2, mapping = aes(xintercept = value)) 
