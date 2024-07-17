"""
Get genotype count from a truth multi-sample VCF
"""
import pysam
import truvari
from collections import Counter
c = Counter()

for entry in pysam.VariantFile("dipcall_bcftools_merge.vcf.gz"):
    if entry.chrom in ['chrX', 'chrY'] or '_' in entry.chrom:
        continue
    for samp in entry.samples.values():
        c[truvari.get_gt(samp['GT']).name] += 1
import json

print(json.dumps(dict(c), indent=4))
