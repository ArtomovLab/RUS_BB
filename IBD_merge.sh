#!/bin/bash

merge="merge-ibd-segments.17Jan20.102.jar"

output='output file'
ibd='file with IBD'
vcf="vcf file of esse + 1000 Genomes"

cat ${ibd} | java -jar ${merge} ${vcf}  plink.GRCh37.map 0.6 1 > ${output}
