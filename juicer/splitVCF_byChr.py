#!/usr/bin/env python

'''
This script will read in a full Sanger mouse genome project SNP VCF file
and separate SNP entries by chromosome into separate VCF files.

Run before snps_from-vcf.sh.

USAGE: python splitVCF_byChr.py  mgp.v2.snps.annot.reformat.vcf

'''

import sys

def processVCF(filepath):
        vcfDict = {}
        vcfDict['header'] = []
        with open(filepath, "r") as f:
                for line in f:
                        if line.startswith('"#') or line.startswith('#'):
                                vcfDict['header'].append(line)
                        else:
                             	lineAry = line.split('\t')
                                chrID = lineAry[0]
                                if chrID not in vcfDict:
                                        vcfDict[chrID] = []
                                vcfDict[chrID].append(line)
        for key in vcfDict:
                output = 'mgp.v2.snps.annot.reformat.chr' + key + '.vcf'
                with open(output, 'w') as g:
                        for item in vcfDict['header']:
                                g.write(item)
                        for item in vcfDict[key]: # for each line in the list
                                g.write(item) # write it to the file

vcfFile = sys.argv[1]
processVCF(vcfFile)
