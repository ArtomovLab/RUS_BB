#!/bin/bash

CHR=$1

beagle="beagle.r1399.jar"

OUTPUT="IBD_chr${CHR}"
vcf="vcf of pruned esse + 1000 Genomes file"

java -Xmx30G -jar ${beagle} ibd=true impute=false gt=${vcf} chrom=${CHR} ibdlod=3 map=plink.GRCh37.map ibdtrim=40 out=${OUTPUT}
