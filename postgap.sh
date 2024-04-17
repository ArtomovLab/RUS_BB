#!/bin/bash

POSTGAP="POSTGAP.py"
database_dir="postgap_libs/databases_dir"

input='input folder'

##RSID

input1='txt file list of rsids'

rsid=${SGE_TASK_ID}
RSID=$(sed "${rsid}q;d" ${input1})
if [[ ! -f ${input}/rsids/SMNE_${RSID}_postgap.txt ]]
then
python ${POSTGAP} --rsID ${RSID} --database_dir ${database_dir} --population EUR --debug --output ${input}/rsids/SMNE_${RSID}_postgap.txt
fi
