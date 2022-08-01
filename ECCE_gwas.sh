#!/bin/bash

#$ -l h_vmem=16G
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=12:00:00
#$ -o gwas.out
#$ -e gwas.err
#$ -cwd
#$ -N gwas
#$ -t 1-2
#$ -tc 2

source /broad/software/scripts/useuse

#use Python-3.6
use Anaconda3
use Java-1.8
use OpenBLAS
use .openblas_for_hail-0.2.20

export PYSPARK_SUBMIT_ARGS="--driver-memory 8G pyspark-shell --executor-memory 8G"
export LD_PRELOAD=/broad/software/free/Linux/redhat_7_x86_64/pkgs/openblas_0.2.20/lib/libopenblas.so

#source /humgen/atgu1/methods/dusoltsev/programms/hail_3_8/bin/activate

pheno=${SGE_TASK_ID}
list=/humgen/atgu1/methods/dusoltsev/biobank/HRC/pheno/pheno.txt
python ECCE_gwas.py -ph ${pheno} -list ${list} 
