"""
How good is canberra similarity at estimating sequence similarity?

I'm going to extract insertions from a VCF and if they're within 10bp, compare all vs all with sequence similarity and
with canberra similarity
"""
import sys
import truvari
import pysam
import itertools
import numpy as np

KMER=4
MAXDIST=500
SIZEMIN=50

def encode_nuc(nuc):
    if nuc == 'A':
        return 0
    if nuc == 'G':
        return 1
    if nuc == 'C':
        return 2
    if nuc == 'G':
        return 3
    return 0

def generate_kmers(sequence, kmer=4):
    result = 0
    for pos, i in enumerate(sequence[:kmer]):
        result += encode_nuc(i) << ((kmer - pos - 1) * 2)
    yield result
    mask = (4**(kmer-1)) - 1
    for i in sequence[1:]:
        nuc = encode_nuc(i)
        result = ((result & mask) << 2) + nuc
        yield result

def kfeat(seq, kmer_len=4):
    """
    Make the kmer array of all kmers and those over min_freq
    """
    ret = np.zeros(4**kmer_len, dtype=np.float32)
    for i in generate_kmers(seq.upper(), kmer_len):
        ret[i] += 1
    return ret


def canberra(a, b, mink=1, maxk=-1):
    """
    canberra similarity of two vectors
    """
    deno = 0.0
    neum = 0.0

    for x, y in zip(a, b):
        total_d = abs(x) + abs(y)
        if total_d > mink:# and (maxk == -1 or total_d < maxk):
            deno += total_d
            neum += abs(x - y)

    if deno == 0.0:
        return 0.0

    if neum == 0.0:
        return 1.0

    return 1.0 - (neum / deno)

def weigh_canberra(a, b, length, mink=1, maxk=-1):
    """
    canberra similarity of two vectors
    """
    deno = 0.0
    neum = 0.0
    def correction(v, corr):
        t = abs(v)
        # Bin by some percent of the length
        #corr = maxk#m_len * maxk
        t = (t / corr) + (t % corr)
        if v < 0:
            return -t
        return t
    view = np.abs(np.hstack([a, b]))
    view = view[view > mink]
    corr = np.median(view)
    thresh = np.percentile(view, 10)
    for x, y in zip(a, b):
        total_d = abs(x) + abs(y)
        if total_d > thresh:
            x = correction(x, corr)
            y = correction(y, corr)
            total_d = abs(x) + abs(y)
        if total_d > mink:
            neum +=  abs(x - y)
            deno += total_d

    if deno == 0.0:
        return 0.0

    if neum == 0.0:
        return 1.0

    return 1.0 - (neum / deno)

def cansim(seq1, seq2, length, kmer=4, mink=1, maxk=-1):
    k1 = kfeat(seq1, kmer)
    k2 = kfeat(seq2, kmer)
    return round(canberra(k1, k2, mink, maxk), 4), round(weigh_canberra(k1, k2, length, mink, maxk), 4)

def analyze_group(cur_group, kmer):
    if len(cur_group) < 2:
        return

    for i, j in itertools.combinations(cur_group, 2):
        if max(truvari.entry_distance(i, j)) > MAXDIST:
            continue
        seq_sim = round(truvari.entry_seq_similarity(i, j), 4)

        can_sim, weigh_can_sim = cansim(i.alts[0], j.alts[0], kmer)

        sz1 = truvari.entry_size(i)
        sz2 = truvari.entry_size(j)
        szbin = truvari.get_sizebin(max(sz1, sz2))
        szsim = round(truvari.sizesim(sz1, sz2)[0], 4)
        #if abs(seq_sim - can_sim) > .5:
            #sys.stderr.write(f"{str(i)}\n{str(j)}\n")
            #print(szbin, szsim, seq_sim, can_sim, sep='\t', file=sys.stderr)
        print(szbin, szsim, seq_sim, can_sim, weigh_can_sim, i.pos == j.pos, sz1, sz2, sep='\t')

if __name__ == '__main__':
    KMER = int(sys.argv[1])
    fn = "../estab_benchmark/giab_comp/GRCh38_HG2-T2TQ100-V1.1.vcf.gz"
    bed_fn = "../estab_benchmark/giab_comp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"

    vcf_i = pysam.VariantFile(fn)
    #bed = truvari.build_region_tree(vcfA=vcf, includebed=bed_fn)
    #vcf_i = truvari.region_filter(vcf, bed)

    print("szbin\tszsim\tseqsim\tcansim\tweigh_can_sim\tsame_start\tsz1\tsz2")

    last_pos = 0
    last_chrom = None
    cur_group = []
    for entry in vcf_i:
        if entry.alts is None or entry.alts[0] == '*':
            continue
        if truvari.entry_variant_type(entry) != truvari.SV.INS:
            continue
        if truvari.entry_size(entry) < SIZEMIN:
            continue
        if last_chrom is None:
            last_chrom = entry.chrom
        if entry.pos - last_pos > MAXDIST or entry.chrom != last_chrom:
            analyze_group(cur_group, KMER)
            cur_group = []
        last_pos = entry.pos
        last_chrom = entry.chrom
        cur_group.append(entry)

    analyze_group(cur_group, KMER)
