#!/bin/bash

#$ -l h_vmem=4G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=300:00:00
#$ -o postgap.out
#$ -e postgap.err
#$ -cwd
#$ -N postgap
#$ -t 1:728290
#$ -tc 100

source /broad/software/scripts/useuse
##use Anaconda3
##use Java-1.8
##pheno="20116_0.gwas.imputed_v3.both_sexes"

##python postgap_preparation.py -pheno ${pheno}

##unuse Anaconda3
##unuse Java-1.8


POSTGAP="/humgen/atgu1/methods/dusoltsev/programms/postgap/POSTGAP.py"
##SUM="/humgen/atgu1/methods/dusoltsev/biobank/ukbb/${pheno}_FILTERED.tsv"
##OUT="/humgen/atgu1/methods/dusoltsev/biobank/ukbb/20116_0.gwas.imputed_v3.both_sexes_POSTGAP1.tsv"
database_dir="/humgen/atgu1/methods/dusoltsev/programms/postgap_libs/databases_dir"

export PYTHONPATH=/humgen/atgu1/methods/dusoltsev/programms/postgap/lib:$PYTHONPATH
export PATH=/humgen/atgu1/methods/dusoltsev/programms/postgap_libs/bin:$PATH
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/humgen/atgu1/methods/dusoltsev/programms/postgap_libs/htslib:/humgen/atgu1/methods/dusoltsev/programms/postgap_libs/libBigWig:/humgen/atgu1/methods/dusoltsev/programms/postgap_libs/gsl/.libs

use Python-2.7
use GCC-5.2
use .gsl-2.6
source /humgen/atgu1/methods/dusoltsev/programms/postgap_libs/postgap2_7/bin/activate
input='/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE'

##python ${POSTGAP} --summary_stats /humgen/atgu1/methods/dusoltsev/biobank/ukbb/${pheno}/${pheno}_chr${CHROM}_FILTERED.tsv --database_dir ${database_dir} --population EUR --debug --output /humgen/atgu1/methods/dusoltsev/biobank/ukbb/${pheno}/${pheno}_chr${CHROM}_POSTGAP.tsv

##CHROM

##CHROM=${SGE_TASK_ID}

##if [[ ${CHROM} == 23 ]]
##then
##CHROM='X'
##elif [[ ${CHROM} == 24 ]]
##then
##CHROM='Y'
##elif [[ ${CHROM} == 25 ]]
##then
##CHROM='MT'
##fi

##python ${POSTGAP} --summary_stats ${input}/SmokingCessation_filtered_chr${CHROM}.txt  --database_dir ${database_dir} --population EUR --debug --output ${input}/SmokingCessation_filtered_chr${CHROM}_postgap.txt 

##test
##python ${POSTGAP} --disease autism --population EUR --database_dir ${database_dir} --output /humgen/atgu1/methods/dusoltsev/biobank/autism.txt

##RSID

input1='/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SMNE/SMNE_filtered_rsids.txt'
rsid=${SGE_TASK_ID}
RSID=$(sed "${rsid}q;d" ${input1})
if [[ ! -f ${input}/rsids/SMNE_${RSID}_postgap.txt ]]
then
python ${POSTGAP} --rsID ${RSID} --database_dir ${database_dir} --population EUR --debug --output ${input}/rsids/SMNE_${RSID}_postgap.txt
fi
