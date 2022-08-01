#! /bin/bash

#$ -cwd
#$ -N fst
#$ -l h_vmem=16G
#$ -q broad
#$ -t 1-29
#$ -tc 29
#$ -l h_rt=24:00:00
#$ -o fst.out
#$ -e fst.err
#$ -pe smp 1
#$ -binding linear:1

source /broad/software/scripts/useuse
use Tabix
use VCFtools

POPs="/humgen/atgu1/methods/dusoltsev/biobank/HRC/fst/fst_ALL.txt"
TARGET="/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS.vcf.bgz"

OUT='regions_WGS_pruned'

id=${SGE_TASK_ID}

POP=$(sed "${id}q;d" ${POPs})

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_SAMARA.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove /humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out1.id
--out ${OUT}/SAMARA_${POP}

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_SPB.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove /humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out1.id
--out ${OUT}/SPB_${POP}

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_ORENBURG.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove /humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out1.id
--out ${OUT}/ORENBURG_${POP}

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_RUS.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove /humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out1.id
--out ${OUT}/RUS_${POP}

N=6
for ((i=1; i <= ${N}; i++))
do
vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_cl${i}.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--out ${OUT}/cl${i}_${POP} 
--remove /humgen/atgu1/methods/dusoltsev/biobank/HRC/plink2.king.cutoff.out1.id 
done
