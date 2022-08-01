#!/bin/bash

#$ -l h_vmem=16G
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=48:00:00
#$ -o qc_plink1.out
#$ -e qc_plink1.err
#$ -cwd
#$ -N qc1

source /broad/software/scripts/useuse

use Anaconda3
use Java-1.8
use R-4.0

export PYSPARK_SUBMIT_ARGS="--driver-memory 8G pyspark-shell --executor-memory 8G"
export LD_PRELOAD=/broad/software/free/Linux/redhat_7_x86_64/pkgs/openblas_0.2.20/lib/libopenblas.so

plink2="/humgen/atgu1/methods/dusoltsev/programms/plink2"

vcf='/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_duprem_all.vcf.gz'
bfile='/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all'
output='/humgen/atgu1/methods/dusoltsev/biobank/HRC'

${plink2} --vcf ${vcf} --king-cutoff 0.08838835 -out all_relatives

${plink2} --bfile ${bfile} --make-king square0

${plink2} -vcf ${vcf} --make-bed --geno 0.01 --hwe 1e-06 --import-dosage-certainty 0.9 --maf 0.01 --out esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS

${plink2} --bfile ${bfile} --indep-pairwise 200 50 0.2 --maf 0.01 --hwe 1e-4 --out ${output}/esse.hrc.dr08_pruned_PLINK2 --exclude ${output}/variants_to_exclude.txt
${plink2} --bfile ${bfile} --extract ${output}/esse.hrc.dr08_pruned_PLINK2.prune.in --out ./esse.hrc.dr08_pruned_PLINK2 --pca 10

rm ./test.txt
for ((i=2;i<=100;i++)); do
Rscript pca_filter.R ${i}
cat ./test.txt >> ./samples_to_remove1.txt
${plink2} --bfile ${bfile} --indep-pairwise 200 50 0.2 --maf 0.01 --hwe 1e-4 --out /humgen/atgu1/methods/dusoltsev/biobank/HRC/pca/esse.hrc.dr08_pruned_PLINK2_${i} --remove ${output}/samples_to_remove1.txt --exclude ${output}/variants_to_exclude.txt
${plink2} --bfile ${bfile} --extract ${output}/pca/esse.hrc.dr08_pruned_PLINK2_${i}.prune.in --remove ${output}/samples_to_remove1.txt --out ./pca/esse.hrc.dr08_pruned_PLINK2_${i} --pca 10 
done

${plink2} --bfile ${bfile} --indep-pairwise 200 50 0.2 --maf 0.01 --hwe 1e-4 --out ${output}/esse.hrc.dr08_pruned_PLINK2_final --remove ${output}/samples_to_remove2.txt --exclude ${output}/variants_to_exclude.txt
${plink2} --bfile ${bfile} --extract ${output}/esse.hrc.dr08_pruned_PLINK2_final.prune.in --out ./esse.hrc.dr08_pruned_PLINK2_final --remove ${output}/samples_to_remove2.txt --pca 10

