#!/bin/bash

#zcat ./ESSE.HRC.vep_annotated_37_CSQ_severe.gz | awk 'BEGIN{FS=OFS="\t"} {print $1":"$2":"$3":"$4,$9}' > ./ESSE.HRC.vep_annotated_37_CSQ_severe_persed.txt
#sed -i '1d' ./ESSE.HRC.vep_annotated_37_CSQ_severe_persed.txt
#sed -i -e '1 s/^/Uploaded_variation\tConsequence\n/;' ./ESSE.HRC.vep_annotated_37_CSQ_severe_persed.txt

zcat ./generated-by-pheweb/sites/sites.tsv | awk 'BEGIN{FS=OFS="\t"} {print $1":"$2":"$3":"$4,$5,$6}' > ./sites1.tsv
sed -i '1d' ./sites1.tsv
sed -i -e '1 s/^/Uploaded_variation\trsid\tnearest_gene\n/;' ./sites1.tsv
