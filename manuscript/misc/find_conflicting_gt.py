import sys
import pysam
import truvari


matcher = truvari.Matcher()
matcher.params.sizefilt = 50
matcher.params.sizemin = 50
matcher.params.sizemax = 10000
matcher.params.chunksize = 100

v = pysam.VariantFile(sys.argv[1])
field = 'DP'

chunks = truvari.chunker(matcher, ('file', v))

n_calls = 0
n_conflict = 0
n_present = 0
for chunk, _ in chunks:
    # Are there multiple deletions
    dels = [_ for _ in chunk['file'] if truvari.entry_variant_type(_) == truvari.SV.DEL and _.alts[0] != '*']
    n_calls += len(dels)
    if len(dels) > 1:
        for i in range(len(dels)-1):
            call1 = dels[i]
            if 1 in call1.samples[0]['GT']:
                n_present += 1
            any_conflict = False
            for j in range(i+1, len(dels)):
                call2 = dels[j]
                ovl = truvari.entry_reciprocal_overlap(call1, call2)
                if ovl and call1.samples[0]['GT'] == (1,1) and call2.samples[0]['GT'] == (1,1):
                    any_conflict = True
            if any_conflict:
                n_conflict += 1
    elif dels:
        if 1 in dels[0].samples[0]['GT']:
            n_present += 1
print(n_calls, n_present, n_conflict, n_conflict / n_present)
