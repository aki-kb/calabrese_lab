#!/bin/bash

#SBATCH -J juicer_diploid6
#SBATCH --mem=2g
#SBATCH --time=48:00:00
#SBATCH -o juicer_diploid6.%j.out
#SBATCH -e juicer_diploid6.%j.err

# Notes: only processes read pairs in which at least on read end overlaps a SNP and the min MAPQ is 10.

# USAGE: sh juicer_diploid_v6.sh
        # REQUIRED input files:
          # CASTEiJ_${i}_paternal_maternal.txt and CASTEiJ_${i}_chr_pos.txt files, outputs from snps_from-vcf.sh
          # merged_nodups.txt file
          # diploid.pl (original script from Juicer)

module load perl

for i in chrX chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
do
      awk -v chr=$i 'OFS="\t" {if (NR > 0) $1="chr"$1; print}' CASTEiJ_${i}_paternal_maternal.txt > CASTEiJ_${i}_paternal_maternal_chr.txt
      awk -v chr=$i '$2==chr && $6==chr && $9 >= 10 && $12 >= 10' merged_nodups.txt | perl diploid.pl -s CASTEiJ_${i}_chr_pos.txt -o CASTEiJ_${i}_paternal_maternal_chr.txt > diploid_${i}.txt
done
