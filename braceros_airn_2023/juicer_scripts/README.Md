
# Juicer-based Hi-C Analysis Pipeline

## Juicer 
```
$ juicer.sh -y mm9_Arima.txt -C 90000000 -L 15840 -Q 15840
```
juicer.sh downloaded from https://github.com/aidenlab/juicer

## Mega 
```
$ mega.sh -g mm9 -L 15840 -Q 15840 <p />
```
mega.sh downloaded from https://github.com/aidenlab/juicer

## Generate chromosome-specific SNP input files for Juicer's diploid.pl
```
$ python splitVCF_byChr.py  mgp.v2.snps.annot.reformat.vcf 
$ sh snps_from-vcf.sh
```
vcftotxt_v6.awk in snps_from-vcf.sh was downloaded and modified from https://github.com/aidenlab/juicer

## Run diploid.pl to generate allelic Hi-C contacts and .hic files
```
$ sh juicer_diploid_v6.sh
$ sh juicer_allelic-make.sh
```
Diploid.pl and commands in juicer_diploid_v6.sh were downloaded and modified from https://github.com/aidenlab/juicer
