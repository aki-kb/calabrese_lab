
# ChIP-, CHART-, and RNA-Seq Analysis

### Generate wiggle file from total, non-allelic mm9-aligned reads
$ perl bigbowtie_to_wig3_mm9.pl B6.sam output color
  
### Extract B6/CAST SNP-overlapping reads
$ perl intersect_reads_snps17.pl B6.sam CAST.sam sanger_mm9_09_14_11 y/n output <p />
sanger_mm9_09_14_11 contains genomic positions of B6/CAST SNPs from Sanger institue <br />
y/n = paired end data?
  
### 40kb binning (w/ 4kb slide) of total, non-allelic read counts
$ bedtools coverage -counts -sorted -g chr_sizes_sort.txt -a mm9_40bin_4slide.bed -b sortedB6.bam output
  
### 10kb binning of allelic read counts
$ perl ase_analyzer8_hDbed.pl input all-chr_mm9_10kb-bin.bed output
