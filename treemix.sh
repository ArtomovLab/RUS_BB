#!/bin/bash

#$ -l h_vmem=16G
#$ -pe smp 2
#$ -binding linear:2
#$ -l h_rt=48:00:00
#$ -o treemix.out
#$ -e treemix.err
#$ -cwd
#$ -N treemix

source /broad/software/scripts/useuse

#use Python-3.6
use Anaconda3

source activate treemix
input='/humgen/atgu1/methods/dusoltsev/biobank/HRC/treemix/ECCE_1000G_all.tsv.gz'

treemix -i ${input} -root YRI -k 500 -bootstrap -o ECCE_1000G
