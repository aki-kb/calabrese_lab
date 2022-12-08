#!/bin/bash

# USAGE: sh snps_from-vcf.sh
	# REQUIRED input files:
		# vftotxt_v6.awk (modified vcftotxt.awk script)
		# mgp.v2.snps.annot.reformat.vcf separated by chr (output from splitVCF_byChr.py)

for i in {1..19} X
do
	awk -f ../vcftotxt_versions/vcftotxt_v6.awk mgp.v2.snps.annot.reformat.chr${i}.vcf
done
