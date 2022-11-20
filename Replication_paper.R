library(ieugwasr)
res <- data.frame(matrix(ncol=12,nrow=0))
colnames(res) <- c('n','chr','position','se' ,'beta', 'p', 'id','rsid','ea','nea','eaf','trait')

d_ff_maf_ <- d_ff_maf %>% filter(AF_RUS > 0.01 & AF_RUS < 0.99)
d_ff_maf_ <- d_ff_maf_  %>% separate('locus',c('CHROM','POS','REF','ALT'),sep='_')
d_ff_maf_ <- d_ff_maf_  %>% unite('locus',c('CHROM','POS'),sep=':')

for (i in  seq(from = 1,to = length(unique(d_ff_maf_$locus)), by = 100)) {
  d <- phewas(
    unique(d_ff_maf_$locus)[i:(i+100)],
    pval = 1e-6,
    batch = c('finn-b'),
    access_token = check_access_token()
  )
  d <- as.data.frame(d)
  print(d)
  print(c(unique(fin$locus)[i],i))
  res <- rbind(res,d)
}
write.table(res,"/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phewas_NEW.txt", sep='\t',quote = F,row.names = F)

fin_phewas <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phewas.txt",sep = '\t', header = T,quote="")
write.table(unique(fin_phewas$trait),"/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phenos_phewas.txt", sep='\t',quote = F,row.names = F)
length(unique(fin_phewas$trait))
length(unique(fin_phewas$rsid))
fin_phewas$trait <- gsub('"','',fin_phewas$trait)


qqman::manhattan(fin_phewas,chr = "chr", bp = "position", p = "p", snp = "rsid",
                 suggestiveline = -log10(5e-8), genomewideline =-log10(5e-8),
                 annotatePval = -log10(5e-8), annotateTop = T,
                 cex.axis=0.8,highlight = potentially_replicated)



fin_phewas <- fin_phewas %>% filter(trait %in% c('Seropositive rheumatoid arthritis, wide',
                                                 'Rheumatoid arthritis (M13_RHEUMA)',
                                                 'Rheumatoid arthritis (M13_RHEUMA_INCLAVO)',
                                                 'Seropositive rheumatoid arthritis, strict definition',
                                                 'Other/unspecified seropositiverheumatoid arthritis',
                                                 'Seropositive rheumatoid arthritis',
                                                 #'Ankylosing spondylitis',
                                                 #'Spondyloarthritis',
                                                 #'Ankylosing spondylitis, strict definition',
                                                 'Hypertension, essential',
                                                 'Hypertension, essential (no controls excluded)',
                                                 'Hypertension',
                                                 'Hypertension (no controls excluded)',
                                                 'Hypertensive diseases (excluding secondary)',
                                                 'Hypertensive diseases',
                                                 'Antihypertensive medication - note that there are other indications',
                                                 'Hypothyroidism,  other/unspecified',
                                                 'Disorders of the thyroid gland',
                                                 'Hypothyroidism, strict autoimmune',
                                                 'Hypothyroidism, strict autoimmune, 3 medication purchases required',
                                                 'Hypothyroidism (congenital or acquired)',
                                                 'Hypothyroidism and >3 levothyroxin purchases',
                                                 'Hypothyroidism, levothyroxin purchases',
                                                 'Hypothyroidism, drug reimbursement',
                                                 'Alzheimer<c3><95>s disease, wide definition',
                                                 'Alzheimer<c3><95>s disease, wide definition (more controls excluded)',
                                                 'Alzheimer disease',
                                                 'Alzheimer disease, including avohilmo',
                                                 'Alzheimer disease (more controls excluded)',
                                                 'Alzheimer<c3><95>s disease (Late onset)',
                                                 'Alzheimer<c3><95>s disease (Late onset) (more controls excluded)',
                                                 'Statin medication',
                                                 'Asthma, including avohilmo',
                                                 'Asthma/COPD (KELA code 203)',
                                                 'Asthma',
                                                 'Asthma (only as main-diagnosis)',
                                                 'Asthma (more controls excluded)',
                                                 'Asthma, unspecified (mode)',
                                                 'Asthma, unspecified (mode) (more controls excluded)',
                                                 'Allergic asthma (mode)',
                                                 'Allergic asthma (mode) (more controls excluded)',
                                                 'Asthma (only as main-diagnosis) (more controls excluded)',
                                                 'Asthma, hospital admissions , main diagnosis only',
                                                 'Sleep apnoea',
                                                 'Sleep apnoea, including avohilmo',
                                                 
                                                 'Type 1 diabetes',
                                                 'Type 1 diabetes, wide definition',
                                                 'Type 1 diabetes, wide definition, subgroup 1',
                                                 'Type 1 diabetes, strict definition',
                                                 'Type 1 diabetes with ophthalmic complications',
                                                 'Type 1 diabetes, strict (exclude DM2)',
                                                 'Type 1 diabetes, strict definition, subgroup 1',
                                                 'Type1 diabetes, definitions combined',
                                                 'Type 1 diabetes with other specified/multiple/unspecified complications',
                                                 'Type 1 diabetes without complications',
                                                 'Type 2 diabetes',
                                                 'Type 2 diabetes, definitions combined, including avohilmo',
                                                 'Diabetes, varying definitions',
                                                 'Diabetes, insuline treatment (Kela reimbursement)',
                                                 'Diabetes, insuline treatment (Kela reimbursement) (more controls excluded)',
                                                 'Other diabetes, wide definition',
                                                 'Type 2 diabetes without complications',
                                                 'Type 2 diabetes with other specified/multiple/unspecified complications',
                                                 'Type 2 diabetes, definitions combined',
                                                 'Type 2 diabetes, strict (exclude DM1)',
                                                 'Type 2 diabetes, wide definition'))

potentially_replicated <- unique(fin_phewas$rsid)
fin_phewas <-  fin_phewas %>% mutate(id = case_when(trait %in% c( 'Diabetes mellitus',
                                                                  'Type 1 diabetes',
                                                                  'Type 1 diabetes, wide definition',
                                                                  'Type 1 diabetes, wide definition, subgroup 1',
                                                                  'Type 1 diabetes, strict definition',
                                                                  'Type 1 diabetes with ophthalmic complications',
                                                                  'Type 1 diabetes, strict (exclude DM2)',
                                                                  'Type 1 diabetes, strict definition, subgroup 1',
                                                                  'Type1 diabetes, definitions combined',
                                                                  'Type 1 diabetes with other specified/multiple/unspecified complications',
                                                                  'Type 1 diabetes without complications',
                                                                  'Type 2 diabetes',
                                                                  'Type 2 diabetes, definitions combined, including avohilmo',
                                                                  'Diabetes, varying definitions',
                                                                  'Diabetes, insuline treatment (Kela reimbursement)',
                                                                  'Diabetes, insuline treatment (Kela reimbursement) (more controls excluded)',
                                                                  'Other diabetes, wide definition',
                                                                  'Type 2 diabetes without complications',
                                                                  'Type 2 diabetes with other specified/multiple/unspecified complications',
                                                                  'Type 2 diabetes, definitions combined',
                                                                  'Type 2 diabetes, strict (exclude DM1)',
                                                                  'Type 2 diabetes, wide definition') ~ 'Diabetes',
                                                    trait %in% c('Seropositive rheumatoid arthritis, wide',
                                                                 'Rheumatoid arthritis (M13_RHEUMA)',
                                                                 'Rheumatoid arthritis (M13_RHEUMA_INCLAVO)',
                                                                 'Seropositive rheumatoid arthritis, strict definition',
                                                                 'Other/unspecified seropositiverheumatoid arthritis',
                                                                 'Seropositive rheumatoid arthritis') ~ 'Arthritis',
                                                    #trait %in% c('Spondyloarthritis',
                                                    #             'Ankylosing spondylitis',
                                                    # 'Ankylosing spondylitis, strict definition') ~ 'Spondylitis',
                                                    trait %in% c('Hypertension, essential',
                                                                 'Hypertension, essential (no controls excluded)',
                                                                 'Hypertension',
                                                                 'Hypertension (no controls excluded)',
                                                                 'Hypertensive diseases (excluding secondary)',
                                                                 'Hypertensive diseases',
                                                                 'Antihypertensive medication - note that there are other indications') ~ 'Hypertension',
                                                    trait %in% c( 'Hypothyroidism,  other/unspecified',
                                                                  'Disorders of the thyroid gland',
                                                                  'Hypothyroidism, strict autoimmune',
                                                                  'Hypothyroidism, strict autoimmune, 3 medication purchases required',
                                                                  'Hypothyroidism (congenital or acquired)',
                                                                  'Hypothyroidism and >3 levothyroxin purchases',
                                                                  'Hypothyroidism, levothyroxin purchases',
                                                                  'Hypothyroidism, drug reimbursement') ~ 'Hypothyroidism',
                                                    trait %in% c( 'Alzheimer<c3><95>s disease, wide definition',
                                                                  'Alzheimer<c3><95>s disease, wide definition (more controls excluded)',
                                                                  'Alzheimer disease',
                                                                  'Alzheimer disease, including avohilmo',
                                                                  'Alzheimer disease (more controls excluded)',
                                                                  'Alzheimer<c3><95>s disease (Late onset)',
                                                                  'Alzheimer<c3><95>s disease (Late onset) (more controls excluded)') ~ 'Alzheimer',
                                                    trait %in% c('Statin medication') ~ 'Statin',
                                                    trait %in% c(  'Sleep apnoea',
                                                                   'Sleep apnoea, including avohilmo') ~ 'Sleep_apnoea',
                                                    trait %in% c('Asthma, including avohilmo',
                                                                 'Asthma/COPD (KELA code 203)',
                                                                 'Asthma',
                                                                 'Asthma (only as main-diagnosis)',
                                                                 'Asthma (more controls excluded)',
                                                                 'Asthma, unspecified (mode)',
                                                                 'Asthma, unspecified (mode) (more controls excluded)',
                                                                 'Allergic asthma (mode)',
                                                                 'Allergic asthma (mode) (more controls excluded)',
                                                                 'Asthma (only as main-diagnosis) (more controls excluded)',
                                                                 'Asthma, hospital admissions , main diagnosis only') ~ 'Asthma'))
all_rus <- unique(fin_phewas$rsid)

for (i in unique(fin_phewas$id)) {
  write.table(fin_phewas %>% filter(id == i),paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phenos_phewas_filtered_",i,".txt",sep=''), sep='\t',quote = F,row.names = F)
}



library('TwoSampleMR')

for (i in unique(fin_phewas$id)) {
  if (i == "Diabetes") {
    #i='Diabetes'
    exposure <- read_exposure_data(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phenos_phewas_filtered_",i,".txt",sep=''),
                                   snp_col='rsid',beta_col='beta',se_col='se',eaf_col='eaf',effect_allele_col='ea',
                                   other_allele_col='nea',	pval_col='p',sep= '\t')
    exposure <- clump_data(exposure,clump_kb = 10000,clump_r2 = 0.1, pop = "EUR")
  }
  else 
  {
    exposure_A <- read_exposure_data(paste("/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phenos_phewas_filtered_",i,".txt",sep=''),
                                     snp_col='rsid',beta_col='beta',se_col='se',eaf_col='eaf',effect_allele_col='ea',
                                     other_allele_col='nea',	pval_col='p',sep= '\t')
    exposure_A <- clump_data(exposure_A,clump_kb = 10000,clump_r2 = 0.1, pop = "EUR")
    exposure <- rbind(exposure,exposure_A)
  }
}

fin_phewas_clumped <- fin_phewas %>% filter(rsid %in% exposure$SNP)
clumped <- unique(fin_phewas$rsid)
fin_phewas_clumped  <- fin_phewas_clumped %>% unite('id_rsid',c('id','rsid'),sep='!')
fin_phewas_clumped  <- fin_phewas_clumped [!duplicated(fin_phewas_clumped$id_rsid),]
fin_phewas_clumped <- fin_phewas_clumped %>% separate('id_rsid',c('id','rsid'),sep='!')
unique(fin_phewas_clumped$rsid)
unique(fin_phewas_clumped$id)

write.table(fin_phewas,"/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phenos_phewas_filtered_clumped.txt", sep='\t',quote = F,row.names = F)

qqman::manhattan(fin_phewas_clumped,chr = "chr", bp = "position", p = "p", snp = "rsid",
                 suggestiveline = -log10(5e-8), genomewideline =-log10(5e-8),
                 annotatePval = -log10(5e-8), annotateTop = T,
                 cex.axis=0.8)

length(table(fin_phewas$id))


#fin_phewas_rus1 <- read.table("/humgen/atgu1/methods/dusoltsev/biobank/HRC/matrix_merge_finnish_enriched_phenos_phewas_filtered_clumped_RUS.txt",sep = '\t', header = T,quote="")

fin_phewas_rus <- fin_phewas_clumped
fin_phewas_rus$p_RUS <- NA
fin_phewas_rus$beta_RUS <- NA
fin_phewas_rus$se_RUS <- NA

fin_phewas_rus[fin_phewas_rus$rsid == 'rs146015702',]$p_RUS <- 8.4e-2 #Sleep_apnoea
fin_phewas_rus[fin_phewas_rus$rsid == 'rs146015702',]$beta_RUS <- 0.59 #Sleep_apnoea
fin_phewas_rus[fin_phewas_rus$rsid == 'rs146015702',]$se_RUS <- 0.34 #Sleep_apnoea

fin_phewas_rus[fin_phewas_rus$rsid == 'rs74630264',]$p_RUS <- 7.0e-1 #Did the doctor tell you that you have / had bronchial asthma? (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs74630264',]$beta_RUS <- -0.18 #Did the doctor tell you that you have / had bronchial asthma? (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs74630264',]$se_RUS <- 0.46 #Did the doctor tell you that you have / had bronchial asthma? (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs187429064' & fin_phewas_rus$id == 'Statin',]$p_RUS <- 5.0e-1 #Combined indicator dyslipidemia or treatment, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs187429064' & fin_phewas_rus$id == 'Statin',]$beta_RUS <- -0.18 #Combined indicator dyslipidemia or treatment, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs187429064' & fin_phewas_rus$id == 'Statin',]$se_RUS <- 0.26 #Combined indicator dyslipidemia or treatment, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs187429064' & fin_phewas_rus$id == 'Diabetes',]$p_RUS <- 1.4e-1 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs187429064' & fin_phewas_rus$id == 'Diabetes',]$beta_RUS <- -1.1 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs187429064' & fin_phewas_rus$id == 'Diabetes',]$se_RUS <- 0.73 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs182611493' & fin_phewas_rus$id == 'Diabetes',]$p_RUS <- 9.7e-2 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs182611493' & fin_phewas_rus$id == 'Diabetes',]$beta_RUS <- -1.7 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs182611493' & fin_phewas_rus$id == 'Diabetes',]$se_RUS <- 1.0 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs182611493' & fin_phewas_rus$id == 'Statin',]$p_RUS <- 6.9e-1 #		Combined indicator dyslipidemia or treatment, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs182611493' & fin_phewas_rus$id == 'Statin',]$beta_RUS <- -0.11 #		Combined indicator dyslipidemia or treatment, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs182611493' & fin_phewas_rus$id == 'Statin',]$se_RUS <- 0.28 #		Combined indicator dyslipidemia or treatment, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs74800719',]$p_RUS <- 4.4e-4 #		Blp
fin_phewas_rus[fin_phewas_rus$rsid == 'rs74800719',]$beta_RUS <- 0.15 #		Blp
fin_phewas_rus[fin_phewas_rus$rsid == 'rs74800719',]$se_RUS <- 0.042 #	 Blp

fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Hypothyroidism',]$p_RUS <- 2.7e-1 #		TSH
fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Hypothyroidism',]$beta_RUS <- -0.14 #		TSH
fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Hypothyroidism',]$se_RUS <- 0.12#	 TSH

fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Diabetes',]$p_RUS <- 7.6e-1 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Diabetes',]$beta_RUS <- -0.11 #		Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Diabetes',]$se_RUS <- 0.34 #	Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Arthritis',]$p_RUS <- 1.5e-1 #	Did your doctor tell you that you have / had Rheumatoid arthritis? (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Arthritis',]$beta_RUS <- -0.64 #		Did your doctor tell you that you have / had Rheumatoid arthritis? (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs56175143' & fin_phewas_rus$id == 'Arthritis',]$se_RUS <- 0.44 #	Did your doctor tell you that you have / had Rheumatoid arthritis? (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs143439093',]$p_RUS <- 2.2e-2 #	Headaches diagnosis: hypertension, migraine or other
fin_phewas_rus[fin_phewas_rus$rsid == 'rs143439093',]$beta_RUS <- -0.32 #		Headaches diagnosis: hypertension, migraine or other
fin_phewas_rus[fin_phewas_rus$rsid == 'rs143439093',]$se_RUS <- 0.14 # Headaches diagnosis: hypertension, migraine or other

fin_phewas_rus[fin_phewas_rus$rsid == 'rs113445611',]$p_RUS <- 4.8e-1 #	Combined indicator dyslipidemia or treatment, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs113445611',]$beta_RUS <- -0.17 #		Combined indicator dyslipidemia or treatment, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs113445611',]$se_RUS <- 0.24 # Combined indicator dyslipidemia or treatment, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs78868334' & fin_phewas_rus$id == 'Diabetes',]$p_RUS <- 1.7e-2 #		Have your relatives had diabetes? Grandfather / grandmother, aunt / uncle, cousins and sisters, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs78868334' & fin_phewas_rus$id == 'Diabetes',]$beta_RUS <- 0.61 #			Have your relatives had diabetes? Grandfather / grandmother, aunt / uncle, cousins and sisters, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs78868334' & fin_phewas_rus$id == 'Diabetes',]$se_RUS <- 0.26 # 	Have your relatives had diabetes? Grandfather / grandmother, aunt / uncle, cousins and sisters, (yes, no)

fin_phewas_rus[fin_phewas_rus$rsid == 'rs77253385',]$p_RUS <- 4.8e-2 #			Combination of hypertension and smoking
fin_phewas_rus[fin_phewas_rus$rsid == 'rs77253385',]$beta_RUS <- 0.16 #				Combination of hypertension and smoking
fin_phewas_rus[fin_phewas_rus$rsid == 'rs77253385',]$se_RUS <- 0.080 # 		Combination of hypertension and smoking

fin_phewas_rus[fin_phewas_rus$rsid == 'rs140000182',]$p_RUS <- 6.2e-1 #				Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs140000182',]$beta_RUS <- -0.18#			Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)
fin_phewas_rus[fin_phewas_rus$rsid == 'rs140000182',]$se_RUS <- 0.37 # 		Self-reported Diabetes mellitus or sugar-lowering therapy, (yes, no)


#fin_phewas_rus <- fin_phewas_rus %>% unite('locus',c('chr','position','nea','ea'))
#fin_phewas_rus <- as.data.frame(merge(fin_phewas_rus,d_ff_maf,by='locus'))
#fin_phewas_rus <- fin_phewas_rus %>% filter(RUS_NEFIN > 0)
qqman::manhattan(fin_phewas_rus,chr = "chr", bp = "position", p = "p_RUS", snp = "rsid",
                 suggestiveline = -log10(0.05/11/8), genomewideline =-log10(0.05/11/8),
                 annotatePval = -log10(0.05/11/8), annotateTop = T,
                 cex.axis=0.8)
