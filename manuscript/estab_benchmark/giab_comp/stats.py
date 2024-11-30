import pysam
import truvari

bed_fn = "HG002_SVs_Tier1_v0.6.bed"
#vcf_fn = "HG002_SVs_Tier1_v0.6.vcf.gz"
vcf_fn = "GIABv0.6_neigh.srt.vcf.gz"

bed_fn = "GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"
vcf_fn = "GIABv1.1_neigh.vcf.gz"

vcf = pysam.VariantFile(vcf_fn)
regions = truvari.build_region_tree(vcf, includebed=bed_fn)

vcf_i = truvari.region_filter(vcf, regions)

# x=has_neigh
# y=het/hom
#check for independence
table = [[0, 0], [0, 0]]
cnt = 0
het = 0
hom = 0
ndel = 0
nins = 0
has_neigh = 0
for entry in vcf_i:
    if truvari.entry_is_filtered(entry):
        continue
    sz = truvari.entry_size(entry)
    if sz < 50 or sz > 10000:
        continue
    if entry.alts[0] in (None, '*'):
        continue
    ty = truvari.entry_variant_type(entry)
    if ty == truvari.SV.DEL:
        ndel += 1
    elif ty == truvari.SV.INS:
        nins += 1
    m_neigh = 0
    if entry.info['NumNeighbors'] != 0:
        has_neigh += 1
        m_neigh = 1
    m_gt = 0
    if entry.samples[0]['GT'] in [(0,1), (1,0)]:
        het += 1
        m_gt = 0
    elif entry.samples[0]['GT'] in [(1,1)]:
        hom += 1
        m_gt = 1
    else:
        continue
    cnt += 1
    table[m_neigh][m_gt] += 1


print(ndel, nins)
print(cnt, het / hom, het, hom, has_neigh, sep='\t')
print(table)
