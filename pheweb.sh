#$ -l h_vmem=16G
#$ -pe smp 1
#$ -binding linear:1
#$ -l h_rt=72:00:00
#$ -o pheweb.out
#$ -e pheweb.err
#$ -cwd
#$ -N pheweb
#$ -t 1:1
#$ -tc 1

source /broad/software/scripts/useuse
export PATH=/humgen/atgu1/methods/dusoltsev/programms/postgap_libs/bin:$PATH
use Anaconda3

python /humgen/atgu1/methods/dusoltsev/biobank/HRC/pheweb_data/pheno_creation.py

input='/humgen/atgu1/methods/dusoltsev/biobank/HRC/gwas_results_all'
output='/humgen/atgu1/methods/dusoltsev/biobank/HRC/pheweb_data/data'
data='/humgen/atgu1/methods/dusoltsev/biobank/HRC/pheno/pheweb.txt'


id=${SGE_TASK_ID}

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
