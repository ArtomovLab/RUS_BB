#!/bin/bash

plink2="plink2"
admixture="admixture"

vcf='vcf file of esse + 1000 Genomes'
bfile='plink bfile of esse + 1000 Genomes'
relatives='txt file with IDs of relatives'

${plink2} -vcf ${vcf} --make-bed --remove ${relatives} --out ${bfile}

${admixture} --cv ${bfile} 8 -j4 --supervised


