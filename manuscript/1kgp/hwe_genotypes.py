import sys
import pysam
import joblib
import numpy as np
import truvari
SZMIN = 50
SZMAX = 10000
def parse_vcf(vcf_file):
    vcf = pysam.VariantFile(vcf_file)
    genotype_counts = []
    for entry in vcf:
        if entry.chrom == 'chrX' or entry.chrom == 'chrY':
            continue
        sz = truvari.entry_size(entry)
        if sz < SZMIN or sz > SZMAX:
            continue

        hom_ref, het, hom_alt, missing = np.uint16(), np.uint16(), np.uint16(), np.uint16()

        for sample in entry.samples.values():
            if None in sample['GT']:
                missing  += 1
                continue
            ac = sample['GT'].count(1)
            if ac == 0:
                hom_ref += 1
            elif ac == 1:
                het += 1
            elif ac == 2:
                hom_alt += 1
        
        genotype_counts.append((hom_ref, het, hom_alt, missing))
    return np.array(genotype_counts, dtype=np.uint16)

if __name__ == '__main__':
    in_vcf = sys.argv[1]
    out_jl = sys.argv[2]
    joblib.dump(parse_vcf(in_vcf), out_jl)
    print('finished', in_vcf)
