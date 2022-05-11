#!/usr/bin/env python

vcfs = "$vcf".split(" ")
with open("vcfList.txt", "w") as fh:
    for vcf in vcfs:
        print(vcf, file=fh)
