
#in python3:
#mt = hl.read_matrix_table('/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc.ht')
#imputed_sex = hl.impute_sex(mt.GT,female_threshold=0.2, male_threshold=0.8)
#imputed_sex.export('/humgen/atgu1/methods/dusoltsev/biobank/HRC/SEX_hail.txt')

library(readxl)
SEX_corrected1 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/SEX_hail.txt",sep= '\t',header =T)

SEX_corrected1 <- SEX_corrected1 %>% mutate(SEX_hail = case_when(is_female == 'true' ~ 2,
                                                                 is_female == 'false' ~ 1))
colnames(SEX_corrected1)[1] <- 'ID'
SEX_corrected1 <- SEX_corrected1 %>% dplyr::select(ID,SEX_hail)
Gen_gwas <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/Phen_gwas.txt",sep = '\t', header = T)
Gen_gwas <- as.data.frame(merge(Gen_gwas,SEX_corrected1,by='ID'))
d <- Gen_gwas %>% dplyr::select(ID,SEX,SEX_hail)
SEX_corrected <- read_excel("/humgen/atgu1/methods/dusoltsev/biobank/ESSE-KBL_manifest_010620_GENDER_FIXED.xlsx",2)
SEX_corrected <- SEX_corrected[!duplicated(SEX_corrected$SampleName),]
colnames(SEX_corrected)[4] <- 'ID'
d <- as.data.frame(merge(d,SEX_corrected,by='ID'))
d <- d %>% mutate(SEX_mismatch = case_when(SEX != SEX_hail ~ 1,
                                            SEX == SEX_hail ~ 0))
table(d$`Sex Mismatch Detected`)
table(d$SEX_mismatch)