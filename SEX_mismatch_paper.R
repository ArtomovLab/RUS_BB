#in python3:

#esse = 'hail mt file of Russian genotypes'
#SEX_hail = 'txt file sex which was difined by hail'
#mt = hl.read_matrix_table(f'{esse}')
#imputed_sex = hl.impute_sex(mt.GT,female_threshold=0.2, male_threshold=0.8)
#imputed_sex.export(f'{SEX_hail}')

library(dplyr)
library(tidyr)
library(ggplot2)
'%!in%' <- function(x,y)!('%in%'(x,y))

SEX_hail <- 'txt file sex which was difined by hail'
esse_pheno <- 'txt file of Russian phenotypes'
GENDER <- 'reference sex'

library(readxl)
SEX_corrected1 <- read.table(SEX_hail,sep= '\t',header =T)

SEX_corrected1 <- SEX_corrected1 %>% mutate(SEX_hail = case_when(is_female == 'true' ~ 2,
                                                                 is_female == 'false' ~ 1))
colnames(SEX_corrected1)[1] <- 'ID'
SEX_corrected1 <- SEX_corrected1 %>% dplyr::select(ID,SEX_hail)
Gen_gwas <- read.table(esse_pheno,sep = '\t', header = T)
Gen_gwas <- as.data.frame(merge(Gen_gwas,SEX_corrected1,by='ID'))
d <- Gen_gwas %>% dplyr::select(ID,SEX,SEX_hail)
SEX_corrected <- read_excel("GENDER",2)
SEX_corrected <- SEX_corrected[!duplicated(SEX_corrected$SampleName),]
colnames(SEX_corrected)[4] <- 'ID'
d <- as.data.frame(merge(d,SEX_corrected,by='ID'))
d <- d %>% mutate(SEX_mismatch = case_when(SEX != SEX_hail ~ 1,
                                            SEX == SEX_hail ~ 0))
table(d$`Sex Mismatch Detected`)
table(d$SEX_mismatch)
