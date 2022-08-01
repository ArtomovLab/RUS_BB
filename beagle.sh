#!/bin/bash

#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=300:00:00
#$ -o refinder.out
#$ -e refinder.err
#$ -cwd
#$ -N refinder
#$ -t 1-22
#$ -tc 22

source /broad/software/scripts/useuse

#use Python-3.6
use Anaconda3
use Java-1.8
use  PLINK2

source /humgen/atgu1/methods/dusoltsev/programms/hail_3_8/bin/activate
CHR=${SGE_TASK_ID}

beagle="/humgen/atgu1/methods/dusoltsev/programms/beagle.r1399.jar"
refinder="/humgen/atgu1/methods/dusoltsev/programms/refined-ibd.17Jan20.102.jar"

OUTPUT="/humgen/atgu1/methods/dusoltsev/biobank/populations/IBD/IBD_1000G_genotype_refinder_chr${CHR}"
vcf="/humgen/atgu1/methods/dusoltsev/biobank/HRC/esse.concat.hrc_qc_duprem_all_pruned_1000G_WGS.vcf.gz"

java -Xmx30G -jar ${beagle} ibd=true impute=false gt=${vcf} chrom=${CHR} ibdlod=3 map=/humgen/atgu1/methods/dusoltsev/biobank/HRC/plink_cm/plink.GRCh37.map ibdtrim=40 out=${OUTPUT}
