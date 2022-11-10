#!/usr/bin/awk -f    
##########
#The MIT License (MIT)
#
# Copyright (c) 2015 Aiden Lab
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#  THE SOFTWARE.
##########
# Parse VCF file into paternal_maternal and chr_pos
# Note that designation is arbitrary as to which comes first
# We'll name "paternal" the first allele and "maternal" the second
#
# Header of VCF
#1. #CHROM
#2. POS
#3. ID
#4. REF
#5. ALT
#6. QUAL
#7. FILTER
#8. INFO
# if GT info present, 9. FORMAT
# 10-end samples
BEGIN {
    OFS="\t";
}
$1 ~ /^#CHROM/ {
    # header line
    if (NF <= 8) {  # NF = number of fields variable
        print "No genotype information available";
        exit;
    }
    else {  # assigning the fields by header line
        for (i=10; i<=NF; i++){
            samples[i]=$i;
        }
    }
    prev="notset";
}
$0 !~ /^#/ {    # $0 = whole input record; telling it to focus only on non-# lines
    # grab the sample fields and parse phasing
    for (i=10; i<=NF; i++){ # start at column 10, end at last column; index each column (i.e. strain)
        split($i, a, ":");  # split(s,a,sep) splits a string 's' into an awk array 'a' using the delimiter 'sep'; split by subfields
        split($4, ref, ""); # split REF
        split($5, alt, ","); # split ALT if more than alternate/non-reference alleles; not necessary for mouse.
        toolong=0;
        for (j=1; j <=length(alt); j++) {
            if (length(alt[j]) > 1) toolong=1;
        }
        # only interested in (unphased) SNPs, not indels
        if (a[1] ~ /\// && length($4)==1 && !toolong) { # if diploid has '/' AND it has more than one alt allele AND it not too long
            # position 0 corresponds to ref SNP
            alt[0]=$4;  # 0th index = REF SNP
            split(a[1], gt, "/");
            # SNP is different
            if (gt[1] == gt[2] && gt[1] == 1 || gt[1] == 2) {   # must be true: a/b, where a=b AND a=1
                if (gt[1]==1) {     # # will take first ALT allele if 1/1
                    print $1":"$2, $4, alt[1] >> samples[i]"_paternal_maternal.txt"; # first allele is "paternal" and listed first
                }
                else {
                    print $1":"$2, $4, alt[2] >> samples[i]"_paternal_maternal.txt";
                }   
                # new chromosome
                if ($1 != prev) {
                    printf("\n%s ", $1) >> samples[i]"_chr_pos.txt";
                prev=$1;
                }
                printf("%d ", $2) >> samples[i]"_chr_pos.txt";
            }
        }
        else {
            print >> samples[i]"_skipped.vcf";
        }
    }
}
END {
    # add newline to file, if it was created
    if (prev != "notset") {
        for (i in samples) {
            printf("\n") >> samples[i]"_chr_pos.txt";
            print samples[i];
        }
    }
}
