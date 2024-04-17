#!/bin/bash

plink2="plink2"

vcf="esse vcf file"
bfile="esse bfile"
output="output directory"
pattern="pattern of filename"

${plink2} --vcf ${vcf} --king-cutoff 0.08838835 -out all_relatives

${plink2} --bfile ${bfile} --make-king square0

${plink2} --bfile ${bfile} --indep-pairwise 200 50 0.2 --maf 0.01 --hwe 1e-4 --out ${output}/${pattern} --exclude ${output}/variants_to_exclude.txt
${plink2} --bfile ${bfile} --extract ${output}/${pattern}.prune.in --out ./${pattern} --pca 10

rm ./test.txt
for ((i=2;i<=100;i++)); do
Rscript pca_filter.R ${i}
cat ./test.txt >> ./samples_to_remove1.txt
${plink2} --bfile ${bfile} --indep-pairwise 200 50 0.2 --maf 0.01 --hwe 1e-4 --out ${pattern}_${i} --remove ${output}/samples_to_remove1.txt --exclude ${output}/variants_to_exclude.txt
${plink2} --bfile ${bfile} --extract ${output}/pca/${pattern}_${i}.prune.in --remove ${output}/samples_to_remove1.txt --out .pca/${pattern}_${i} --pca 10 
done

${plink2} --bfile ${bfile} --indep-pairwise 200 50 0.2 --maf 0.01 --hwe 1e-4 --out ${output}/${pattern}_final --remove ${output}/samples_to_remove2.txt --exclude ${output}/variants_to_exclude.txt
${plink2} --bfile ${bfile} --extract ${output}/${pattern}_final.prune.in --out ./${pattern}_final --remove ${output}/samples_to_remove2.txt --pca 10

