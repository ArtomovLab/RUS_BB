library(magrittr)
inputFile <- "/humgen/atgu1/methods/dusoltsev/biobank/new_ecce/vep/ESSE.HRC.vep_annotated_37.vcf"
outputFile <- "/humgen/atgu1/methods/dusoltsev/biobank/new_ecce/vep/ESSE.HRC.vep_annotated_37_CSQ"
##annotation rating
annotationRanking <-
  c(
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "synonymous_variant",
    "coding_sequence_variant",
    "mature_miRNA_variant",
    "5_prime_UTR_variant",
    "3_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intron_variant",
    "NMD_transcript_variant",
    "non_coding_transcript_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "TFBS_ablation",
    "TFBS_amplification",
    "TF_binding_site_variant",
    "regulatory_region_ablation",
    "regulatory_region_amplification",
    "feature_elongation",
    "regulatory_region_variant",
    "feature_truncation",
    "missense_variant&splice_region_variant",
    "splice_region_variant&synonymous_variant",
    "splice_region_variant&intron_variant",
    "splice_region_variant&intron_variant&NMD_transcript_variant",
    "splice_region_variant&intron_variant&non_coding_transcript_variant",
    "splice_region_variant&non_coding_transcript_exon_variant",
    "intron_variant&NMD_transcript_variant",
    "intron_variant&non_coding_transcript_variant",
    "3_prime_UTR_variant&NMD_transcript_variant",
    "transcribed_unitary_pseudogene",
    "intergenic_variant"
  )

incon <- gzcon(file(inputFile, open="rb"))
outConSevere <- file(paste(outputFile, "_severe", sep = ""), "w")

while ( TRUE ) {

  line = readLines(incon, n = 1)
  line
  if ( length(line) == 0 ) {
    break
  }
  if (grepl(pattern = "##", x = line) ){
    next
  }
  if (grepl(pattern = "INFO", x = line) ){
    line <- strsplit(x = line, split = "\t") %>% unlist
    writeLines(paste('CHROM','POS','REF','ALT','ID','FILTER','Gene','HGNC_ID', 'Consequence','HGVSp', 'LoF', 'LoF_filter', 'LoF_flags', 'LoF_info', 'consequence_field', sep = "\t"),
               con = outConSevere)
    next
  }
  originalLine <- line
  line <- strsplit(x = line, split = "\t") %>% unlist
  CSQ_field <- sub('.*CSQ=', '', line[8]) 
  CSQ_field <- paste('CSQ=',CSQ_field,sep='') 
  CSQ_field1 <- CSQ_field
  CSQ_field <- strsplit(x = CSQ_field, split = ",") %>% unlist
  
  effect <- lapply(CSQ_field, function(x){
    strsplit(x, split = "\\|") %>%
      unlist %>% .[2] }) %>%  unlist
  
  gene <- lapply(CSQ_field, function(x){
    strsplit(x, split = "\\|") %>%
      unlist %>% .[4] }) %>% unlist
  
  id <- lapply(CSQ_field, function(x){
    strsplit(x, split = "\\|") %>%
      unlist %>% .[5] }) %>% unlist
  
  outputPrefix <- strsplit(x = originalLine, split = "\t") %>%
    unlist %>% .[c(1,2,4,5,3,7)] %>% paste(., collapse = "\t")
  print(outputPrefix)
  #for (i in 1:length(effect)){
  #  writeLines(paste(outputPrefix, effect[i], gene[i],, sep = "\t"), 
  #             con = outConAll)
  #}                                                  
  severeImpact <- min(which(annotationRanking %in%
                              effect))
  
  severeImpact <- unique(annotationRanking[severeImpact])
  severeGene <- unique(gene[which(effect == severeImpact)])
  
  severeId <- unique(id[which(effect == severeImpact)])
  #severeAA_pos <- unique(aa_pos[which(effect == severeImpact)])
  #severeAA_cange <- unique(aa[which(effect == severeImpact)])
  
  for(i in 1:length(severeGene)){
    writeLines(paste(outputPrefix,severeGene[i],severeId[i],severeImpact,'','','','','',CSQ_field1, sep = "\t"),
               con = outConSevere)
  }
}

close(incon)
close(outConSevere)


