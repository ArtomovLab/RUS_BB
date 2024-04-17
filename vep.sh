#!/bin/bash

target="input folder"
vep="ensembl-vep/vep"
cache=".vep"
vcf="esse vcf file for annotation"

chr=$1

if [[ ${chr} == 23 ]]
then
chr='X'
fi

tabix -h ${target}/{vcf}.vcf.bgz ${chr} > ${vcf}.chr${chr}.vcf

${vep} --cache --dir_cache ${cache} --cache_version 105 --grch37  -i ${target}/${vcf}.chr${chr}.vcf -o ${vcf}_annotated.chr${chr}.vcf --force_overwrite --vcf

bcftools concat -o ${vcf}_annotated_37.vcf ${vcf}_annotated.chr1.vcf ${vcf}_annotated.chr2.vcf ${vcf}_annotated.chr3.vcf ${vcf}_annotated.chr4.vcf ${vcf}_annotated.chr5.vcf ${vcf}_annotated.chr6.vcf ${vcf}_annotated.chr7.vcf ${vcf}_annotated.chr8.vcf ${vcf}_annotated.chr9.vcf ${vcf}_annotated.chr10.vcf ${vcf}_annotated.chr11.vcf ${vcf}_annotated.chr12.vcf ${vcf}_annotated.chr13.vcf ${vcf}_annotated.chr14.vcf ${vcf}_annotated.chr15.vcf ${vcf}_annotated.chr16.vcf ${vcf}_annotated.chr17.vcf ${vcf}_annotated.chr18.vcf ${vcf}_annotated.chr19.vcf ${vcf}_annotated.chr20.vcf ${vcf}_annotated.chr21.vcf ${vcf}_annotated.chr22.vcf ${vcf}_annotated.chrX.vcf
 
Rscript ~/CSQ.R
