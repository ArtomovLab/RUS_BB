
#!/bin/bash

#$ -l h_vmem=8G
#$ -pe smp 4
#$ -binding linear:4
#$ -l h_rt=300:00:00
#$ -o admixture1.out
#$ -e admixture1.err
#$ -cwd
#$ -N admixture1

source /broad/software/scripts/useuse

#use Python-3.6
use Anaconda3
use Java-1.8
use Tabix
use Bcftools
BCFTOOLS="/humgen/atgu1/methods/nkolosov/tools/bcftools/bcftools"

plink2="/humgen/atgu1/methods/dusoltsev/programms/plink2"
admixture="/humgen/atgu1/methods/dusoltsev/programms/dist/admixture_linux-1.3.0/admixture"

output='/humgen/atgu1/methods/dusoltsev/biobank/HRC'
vcf='/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_annotated.vcf.bgz'
bfile='/humgen/atgu1/methods/dusoltsev/biobank/HRC/ADMIXTURE/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS'

${BCFTOOLS} annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' -Ou ${output}/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS.vcf.bgz -Oz -o ${output}/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS_annotated.vcf.bgz
${plink2} -vcf ${vcf} --make-bed --remove ${output}/all_relatives.king.cutoff.out.id --out ${output}/ADMIXTURE/pop9/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS1

${admixture} --cv ${output}/ADMIXTURE/pop9/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS1.bed 8 -j4 --supervised


