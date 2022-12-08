#!/bin/bash

#SBATCH -J juicer_allelic
#SBATCH --mem=48g
#SBATCH -n 2
#SBATCH --time=48:00:00
#SBATCH -o juicer_allelic.%j.out
#SBATCH -e juicer_allelic.%j.err

# USAGE: sbatch juicer_allelic-make.sh
        # run in directory with juicer_diploid_v6.sh output files.

juiceDir="/juicer/scripts"

## concatenate all chromosome-specific diploid.pl output files (from juicer_diploid_v6.sh).
cat diploid_chr1.txt diploid_chr2.txt diploid_chr3.txt diploid_chr4.txt diploid_chr5.txt diploid_chr6.txt diploid_chr7.txt diploid_chr8.txt diploid_chr9.txt diploid_chr10.txt diploid_chr11.txt diploid_chr12.txt diploid_chr13.txt diploid_chr14.txt diploid_chr15.txt diploid_chr16.txt diploid_chr17.txt diploid_chr18.txt diploid_chr19.txt diploid_chrX.txt > diploid_all-chr.txt

## split allelic hi-c contacts by maternal (CAST-assigned) and paternal (B6-assigned)
#awk -f ${juiceDir}/diploid_split.awk diploid_all-chr.txt

## prepare hic file of maternal (CAST-assigned) hic contacts.
        ## for B/C TSCs (e.g. 10-1), maternal is still CAST, so need to rename to paternal.
sort -k2,2d -m maternal.txt maternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > CASTEiJ_combined.txt
awk '{if ($3 > $7){ print $1, $6, $7, $8, $9, $11, $2, $3, $4, $5, $10}else {print}}' CASTEiJ_combined.txt > CASTEiJ_combined2.txt     ## sort columns so that first read end is less than second read end
sort -k3,3d -k7,7d CASTEiJ_combined2.txt > CASTEiJ_combined_output.txt         ## sort by chr and start
java -Xmx4048m -jar ${juiceDir}/juicer_tools.jar pre CASTEiJ_combined_output.txt CASTEiJ_combined_output.hic mm9

## prepare hic file of paternal (B6-assigned) hic contacts.
        ## for B/C TSCs (e.g. 10-1), paternal is still B6, so need to rename to maternal.
sort -k2,2d -m paternal.txt paternal_both.txt | awk '{split($11,a,"/"); print a[1], $1,$2,$3,$4,$5,$6,$7,$8,"100","100"}' > B6_combined.txt
awk '{if ($3 > $7){ print $1, $6, $7, $8, $9, $11, $2, $3, $4, $5, $10}else {print}}' B6_combined.txt > B6_combined2.txt
sort -k3,3d -k7,7d B6_combined2.txt > B6_combined_output.txt
java -Xmx4048m -jar ${juiceDir}/juicer_tools.jar pre B6_combined_output.txt B6_combined_output.hic mm9
