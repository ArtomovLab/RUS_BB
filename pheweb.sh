#!/bin/bash

input='path to gwas results'
output='output folder'
data='txt file with list of phenotypes'

id=$1

gwas=$(sed "${id}q;d" ${data})

LIN=${input}/${gwas}.txt.bgz
LIN1=${input}/${gwas}_logistic.txt.bgz

if [[ -f $LIN ]]
then

zcat ${input}/${gwas}*.txt.bgz | awk 'BEGIN{FS=OFS="\t"} {gsub(/,/, "!", $2)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/:/, "!", $1)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/,/, "!", $4)} 1'  >  ${output}/${gwas}.phe
awk -F"!" '$1=$1' OFS="\t" ${output}/${gwas}.phe > ${output}/${gwas}.phe1
awk 'BEGIN{FS=OFS="\t"} {if ($1 != "X") {print $1,$2,$3,$4,$5,$7,$8,$11,$12,$14}}' ${output}/${gwas}.phe1 > ${output}/${gwas}.phe
sed -i "s/\[//g" ${output}/${gwas}.phe
sed -i "s/]//g" ${output}/${gwas}.phe
sed -i 's/"//g' ${output}/${gwas}.phe
sed -i '1d' ${output}/${gwas}.phe

sed -i -e '1 s/^/chrom\tpos\tref\talt\trsid\taf\tn\tbeta\tse\tpval\n/;' ${output}/${gwas}.phe
rm ${output}/${gwas}.phe1
bgzip ${output}/${gwas}.phe

fi

if [[ -f $LIN1 ]]
then

zcat ${input}/${gwas}_logistic.txt.bgz | awk 'BEGIN{FS=OFS="\t"} {gsub(/,/, "!", $2)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/:/, "!", $1)} 1' | awk 'BEGIN{FS=OFS="\t"} {gsub(/,/, "!", $4)} 1' >  ${output}/${gwas}_logistic.phe
awk -F"!" '$1=$1' OFS="\t" ${output}/${gwas}_logistic.phe > ${output}/${gwas}_logistic.phe1

awk 'BEGIN{FS=OFS="\t"} {if ($1 != "X") {print $1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12,$14}}' ${output}/${gwas}_logistic.phe1 > ${output}/${gwas}_logistic.phe
sed -i "s/\[//g" ${output}/${gwas}_logistic.phe
sed -i "s/]//g" ${output}/${gwas}_logistic.phe
sed -i 's/"//g' ${output}/${gwas}_logistic.phe
sed -i '1d' ${output}/${gwas}_logistic.phe
sed -i -e '1 s/^/chrom\tpos\tref\talt\trsid\taf\tn\tn_controls\tn_cases\tbeta\tse\tpval\n/;' ${output}/${gwas}_logistic.phe
rm ${output}/${gwas}_logistic.phe1
bgzip  ${output}/${gwas}_logistic.phe

fi
