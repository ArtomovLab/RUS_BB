#!/bin/bash

IBDne="ibdne.23Apr20.ae9.jar"
POPs="txt file with populations"

id=$1

POP=$(sed "${id}q;d" ${POPs})

output='output file'
ibd='IBD file'
vcf="vcf file esse + 1000 Genomes"

mincm=2
cat ${ibd}_${POP}.ibd | java -Xmx30G -jar ${IBDne} mincm=${mincm} map=plink.GRCh37.map out=${output}_${mincm}_${POP}

