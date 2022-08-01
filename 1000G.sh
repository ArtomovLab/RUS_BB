#!/bin/bash

#$ -l h_vmem=16G
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=48:00:00
#$ -o 1000G.out
#$ -e 1000G.err
#$ -cwd
#$ -N G1000

source /broad/software/scripts/useuse

#use Python-3.6
use Anaconda3
use Java-1.8
use OpenBLAS
use .openblas_for_hail-0.2.20

export PYSPARK_SUBMIT_ARGS="--driver-memory 8G pyspark-shell --executor-memory 8G"
export LD_PRELOAD=/broad/software/free/Linux/redhat_7_x86_64/pkgs/openblas_0.2.20/lib/libopenblas.so

#source /humgen/atgu1/methods/dusoltsev/biobank/HRC/hail/bin/activate

python 1000G.py
