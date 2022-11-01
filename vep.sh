#!/bin/bash

#$ -l h_vmem=16G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=72:00:00
#$ -o vep.out
#$ -e vep.err
#$ -cwd
#$ -N vep
##$ -t 1-23
##$ -tc 23

source /broad/software/scripts/useuse
export PATH=/humgen/atgu1/methods/dusoltsev/programms/postgap/bin:$PATH

use MySQL-5.6
use Bcftools
use R-4.0

export PATH=/humgen/atgu1/methods/dusoltsev/programms/perl-5.16.2:$PATH
##export PATH=/humgen/atgu1/methods/dusoltsev/programms/perl-5.16.2/lib:$PATH

target="/humgen/atgu1/methods/dusoltsev/biobank/new_ecce/vep"
vep="/humgen/atgu1/methods/dusoltsev/programms/ensembl-vep/vep"
cache="/humgen/atgu1/methods/dusoltsev/programms/.vep"

chr=${SGE_TASK_ID}

if [[ ${chr} == 23 ]]
then
chr='X'
fi

tabix -h ${target}/ESSE.HRC.vep.vcf.gz ${chr} > ESSE.HRC.vep.chr${chr}.vcf

${vep} --cache --dir_cache ${cache} --cache_version 105 --grch37  -i ${target}/ESSE.HRC.vep.chr${chr}.vcf -o ESSE.HRC.vep_annotated.chr${chr}.vcf --force_overwrite --vcf

bcftools concat -o ESSE.HRC.vep_annotated_37.vcf ESSE.HRC.vep_annotated.chr1.vcf ESSE.HRC.vep_annotated.chr2.vcf ESSE.HRC.vep_annotated.chr3.vcf ESSE.HRC.vep_annotated.chr4.vcf ESSE.HRC.vep_annotated.chr5.vcf ESSE.HRC.vep_annotated.chr6.vcf ESSE.HRC.vep_annotated.chr7.vcf ESSE.HRC.vep_annotated.chr8.vcf ESSE.HRC.vep_annotated.chr9.vcf ESSE.HRC.vep_annotated.chr10.vcf ESSE.HRC.vep_annotated.chr11.vcf ESSE.HRC.vep_annotated.chr12.vcf ESSE.HRC.vep_annotated.chr13.vcf ESSE.HRC.vep_annotated.chr14.vcf ESSE.HRC.vep_annotated.chr15.vcf ESSE.HRC.vep_annotated.chr16.vcf ESSE.HRC.vep_annotated.chr17.vcf ESSE.HRC.vep_annotated.chr18.vcf ESSE.HRC.vep_annotated.chr19.vcf ESSE.HRC.vep_annotated.chr20.vcf ESSE.HRC.vep_annotated.chr21.vcf ESSE.HRC.vep_annotated.chr22.vcf ESSE.HRC.vep_annotated.chrX.vcf
 
Rscript ~/CSQ.R
