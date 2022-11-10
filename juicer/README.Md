
# Hi-C Analysis Pipeline

## Juicer 
$ juicer.sh -y mm9_Arima.txt -C 90000000 -L 15840 -Q 15840

## Mega 
#### $ mega.sh -g mm9 -L 15840 -Q 15840

## Generate chromosome-specific SNP input files for Juicer's diploid.pl
$ module load python/2.7.3 <br />
$ python splitVCF_byChr.py <sanger-mouse-snps>.vcf <br />
$ sh snps_from-vcf.sh

## Generate allelic Hi-C contacts and .hic files
$ sh juicer_diploid_v6.sh <br />
$ sh juicer_allelic-make.sh

## Viewpoint
$ java -Xmx4048m -jar juicer_tools.jar pre merged_nodups.txt output.hic mm9
