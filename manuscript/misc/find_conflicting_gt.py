import sys
import pysam
import truvari


matcher = truvari.Matcher()
matcher.params.sizefilt = 50
matcher.params.sizemin = 50
matcher.params.sizemax = 10000
matcher.params.chunksize = 100

v = pysam.VariantFile(sys.argv[1])

chunks = truvari.chunker(matcher, ('file', v))
def conflict(call1, call2):
    ovl = truvari.entry_reciprocal_overlap(call1, call2, ins_inflate=False)
    if not ovl:
        return False
    samp = 0
    while samp < len(call1.samples):
        c1_gt = truvari.get_gt(call1.samples[samp]['GT'])
        c2_gt = truvari.get_gt(call2.samples[samp]['GT'])
        if c1_gt == truvari.GT.HOM and c2_gt in [truvari.GT.HOM, truvari.GT.HET]:
            #print(ovl, samp)
            return True
        if c2_gt == truvari.GT.HOM and c1_gt in [truvari.GT.HOM, truvari.GT.HET]:
            #print(ovl, samp)
            return True
        samp += 1
    return False

n_calls = 0
n_conflict = 0
for chunk, _ in chunks:
    # Are there multiple deletions
    # This is only checking deletions...
    variants = chunk['file']
    tcnt = len(variants)
    n_calls += tcnt
    if tcnt == 1:
        continue
    status = [False] * len(variants)
    for i in range(tcnt - 1):
        if status[i]:
            continue
        call1 = variants[i]
        for j in range(i+1, tcnt):
            if status[j]:
                continue
            call2 = variants[j]
            if conflict(call1, call2):
                #print(call1, call2)
                status[i] = True
                status[j] = True
    n_conflict += sum(status)
print(n_calls, n_conflict, n_conflict / n_calls)
