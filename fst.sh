#! /bin/bash

POPs="txt file with populations"
TARGET="pruned vcf file esse + 1000 Genomes"
relatives="IDs of esse relatives txt file"

OUT='regions_WGS_pruned'

id=$1

POP=$(sed "${id}q;d" ${POPs})

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_SAMARA.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove ${relatives}
--out ${OUT}/SAMARA_${POP}

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_SPB.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove ${relatives}
--out ${OUT}/SPB_${POP}

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_ORENBURG.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove ${relatives}
--out ${OUT}/ORENBURG_${POP}

vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_RUS.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--remove ${relatives}
--out ${OUT}/RUS_${POP}

N=6
for ((i=1; i <= ${N}; i++))
do
vcftools --gzvcf  $TARGET \
--weir-fst-pop samples_WGS_cl${i}.txt \
--weir-fst-pop samples_WGS_${POP}.txt \
--out ${OUT}/cl${i}_${POP} 
--remove ${relatives} 
done
