
# Hi-C Analysis Pipeline

## Juicer 
$ juicer.sh -y mm9_Arima.txt -C 90000000 -L 15840 -Q 15840 <p />
juicer.sh downloaded from https://github.com/aidenlab/juicer

## Mega 
$ mega.sh -g mm9 -L 15840 -Q 15840 <p />
mega.sh downloaded from https://github.com/aidenlab/juicer

## Generate chromosome-specific SNP input files for Juicer's diploid.pl
$ python splitVCF_byChr.py <sanger-mouse-snps>.vcf <br />
$ sh snps_from-vcf.sh

## Run diploid.pl to generate allelic Hi-C contacts and .hic files
$ sh juicer_diploid_v6.sh <br />
$ sh juicer_allelic-make.sh <p />
Diploid.pl and commands in juicer_diploid_v6.sh were downloaded and modified from https://github.com/aidenlab/juicer
