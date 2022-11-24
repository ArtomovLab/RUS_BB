#!/bin/bash

#$ -l h_vmem=8G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=12:00:00
#$ -o gprior.out
#$ -e gprior.err
#$ -cwd
#$ -N gprior

source /broad/software/scripts/useuse

use Python-3.6

source /humgen/atgu1/methods/dusoltsev/programms/GPrior/vprior/bin/activate


pheno='20116_0.gwas.imputed_v3.both_sexes'

##for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X

##do
##GPRIOR_PRE="/humgen/atgu1/methods/dusoltsev/programms/GPrior/process_postgap.py"
##POSTGAP="/humgen/atgu1/methods/dusoltsev/biobank/ukbb/${pheno}/${pheno}_chr${i}_POSTGAP.tsv"
##GPRIOR="/humgen/atgu1/methods/dusoltsev/biobank/ukbb/${pheno}/${pheno}_chr${i}_GPRIOR.tsv"

##python process_postgap.py -i ${POSTGAP} -o ${GPRIOR}
##done

input='/humgen/atgu1/methods/dusoltsev/biobank/POSTGAP/SmokingCessation'
##python gprior.py -i ${input}/${pheno}/${pheno}_GPRIOR.tsv -ts ${input}/${pheno}_GPRIOR_train.tsv -o ${input}/${pheno}_GPRIOR_results.tsv

python process_postgap.py -i ${input}/SmokingCessation_postgap.tsv -o ${input}/SmokingCessation_postgap_gprior.tsv
python gprior.py -i ${input}/SmokingCessation_postgap_gprior.tsv -ts ${input}/SmokingCessation_train.tsv -o ${input}/SmokingCessation_gprior_results.txt
