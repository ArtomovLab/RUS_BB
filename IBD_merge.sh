#!/bin/bash

#$ -l h_vmem=16G
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=24:00:00
#$ -o IBD_merge.out
#$ -e IBD_merge.err
#$ -cwd
#$ -N IBD_merge
##$ -t 1-32
##$ -tc 32

source /broad/software/scripts/useuse

#use Python-3.6
use Anaconda3
use Java-1.8
use Tabix
use Bcftools
use VCFtools
BCFTOOLS="/humgen/atgu1/methods/nkolosov/tools/bcftools/bcftools"

plink2="/humgen/atgu1/methods/dusoltsev/programms/plink2"
IBDne="/humgen/atgu1/methods/dusoltsev/programms/ibdne.23Apr20.ae9.jar"
merge="/humgen/atgu1/methods/dusoltsev/programms/merge-ibd-segments.17Jan20.102.jar"
POPs="/humgen/atgu1/methods/dusoltsev/biobank/HRC/ibdNE/POPs.txt"

##id=${SGE_TASK_ID}

##POP=$(sed "${id}q;d" ${POPs})

POP='FIN'
output='/humgen/atgu1/methods/dusoltsev/biobank/HRC'
ibd='/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/IBD_1000G_genotype'
vcf="/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS.vcf.gz"

cat ${ibd}.ibd | java -jar ${merge} ${vcf}  ${output}/plink_cm/plink.GRCh37.map 0.6 1 > /humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/IBD_1000G_genotype_merge.ibd
