#!/bin/bash

input='input folder'

python process_postgap.py -i ${input}/SMNE_postgap.tsv -o ${input}/SMNE_postgap_gprior.tsv
python gprior.py -i ${input}/SMNE_postgap_gprior.tsv -ts ${input}/SMNE_train.tsv -o ${input}/SMNE_gprior_results.txt
