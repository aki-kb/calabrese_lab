
# Hi-C Analysis Pipeline

## Juicer 
$ juicer.sh -y mm9_Arima.txt -C 90000000 -L 15840 -Q 15840 <p />
juicer.sh downloaded from https://github.com/aidenlab/juicer

## Mega 
$ mega.sh -g mm9 -L 15840 -Q 15840 <p />
mega.sh downloaded from https://github.com/aidenlab/juicer

## Generate chromosome-specific SNP input files for Juicer's diploid.pl
$ module load python/2.7.3 <br />
$ python splitVCF_byChr.py <sanger-mouse-snps>.vcf <br />
$ sh snps_from-vcf.sh

## Generate allelic Hi-C contacts and .hic files
$ sh juicer_diploid_v6.sh <br />
$ sh juicer_allelic-make.sh <p />
juicer_diploid_v6.sh downloaded and modified from https://github.com/aidenlab/juicer
  
## Viewpoint
$ java -Xmx4048m -jar juicer_tools.jar pre merged_nodups.txt output.hic mm9 <p />
juicer_tools.jar (version 1.9.9) downloaded from https://github.com/aidenlab/juicer
