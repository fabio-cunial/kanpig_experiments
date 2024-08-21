"""
Count how many variants have at least one neighbor which is:
* within 500bp
* of the same type
* over 95% size and sequence similar
"""

import sys
import pysam
import truvari
vcf = pysam.VariantFile(sys.argv[1])

matcher = truvari.Matcher()
matcher.params.sizemin = 50
matcher.params.sizefilt = 50
matcher.params.sizemax = 10000
matcher.params.refdist = 500
matcher.params.seqsim = 0.95
matcher.params.sizesim = 0.95
matcher.params.pick = 'multi'

chunks = truvari.chunker(matcher, ('base', vcf))


matched = 0
checked = 0
for m, _ in chunks:
    if not m['base'] or m['base'][0].chrom in ['chrX', 'chrY']:
        continue
    for i, call in enumerate(m['base']):
        checked += 1
        j = 0
        is_matched = False
        while not is_matched and j < len(m['base']):
            if i != j:
                mat = matcher.build_match(call, m['base'][j], skip_gt=True, short_circuit=True)
                is_matched = mat.state
            j += 1
        if is_matched:
            matched += 1

print(checked, matched)
