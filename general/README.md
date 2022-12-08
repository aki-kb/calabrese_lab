
# Custom Perl Scripts <br />
Scripts used for ChIP-/CHART-/RNA-Seq analyses.

### Generate wiggle file from total, non-allelic mm9-aligned reads
```
$ perl bigbowtie_to_wig3_mm9.pl B6.sam output color
```
B6.sam contains MAPQ>=30 reads <p />
  
### Extract B6/CAST SNP-overlapping reads
```
$ perl intersect_reads_snps17.pl B6.sam CAST.sam sanger_mm9_09_14_11 y/n output
```
< br / >
B6.sam and CAST/sam contain MAPQ>=30 reads <br />
sanger_mm9_09_14_11 contains mm9 genomic positions of B6/CAST SNPs from Sanger institue <br />
y/n = paired end data?
  
### 40kb binning (w/ 4kb slide) of total, non-allelic read counts
```
$ bedtools coverage -counts -sorted -g chr_sizes_sort.txt -a mm9_40bin_4slide.bed -b sortedB6.bam output
```
chr_sizes_sort.txt contains mm9 chr sizes <br />
mm9_40bin_4slide.bed contains 40kb-sized bins tiled across each mm9 chromosome every 4kb. <br />
sortedB6.bam contains MAPQ>=30 reads


### 10kb binning of allelic read counts
```
$ perl ase_analyzer8_hDbed.pl input all-chr_mm9_10kb-bin.bed output
```
input is the output file from intersect_reads_snps17.pl <br />
all-chr_mm9_10kb-bin.bed contains 10kb-sized bins tiled across each mm9 chromosome every 10kb.
