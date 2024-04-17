#!/bin/bash

zcat ./generated-by-pheweb/sites/sites.tsv | awk 'BEGIN{FS=OFS="\t"} {print $1":"$2":"$3":"$4,$5,$6}' > ./sites1.tsv
sed -i '1d' ./sites1.tsv
sed -i -e '1 s/^/Uploaded_variation\trsid\tnearest_gene\n/;' ./sites1.tsv
