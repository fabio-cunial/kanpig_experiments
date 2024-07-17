import pysam
import truvari

bed_fn = "HG002_SVs_Tier1_v0.6.bed"
vcf_fn = "HG002_SVs_Tier1_v0.6.vcf.gz"

bed_fn = "/Users/english/code/truvari/tickets/sawfish/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed"
vcf_fn = "/Users/english/code/truvari/tickets/sawfish/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz"

vcf = pysam.VariantFile(vcf_fn)
regions = truvari.build_region_tree(vcf, includebed=bed_fn)

vcf_i = truvari.region_filter(vcf, regions)

cnt = 0
het = 0
hom = 0

for entry in vcf_i:
    if truvari.entry_is_filtered(entry):
        continue
    sz = truvari.entry_size(entry)
    if sz < 50 or sz > 50000:
        continue
    if entry.samples[0]['GT'] in [(0,1), (1,0)]:
        het += 1
    elif entry.samples[0]['GT'] in [(1,1)]:
        hom += 1
    else:
        continue
    cnt += 1
print(cnt, het / hom, het, hom, sep='\t')
