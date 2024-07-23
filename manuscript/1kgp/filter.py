import sys
import pysam
import truvari

def get_exchet(entry):
    if "ExcHet" not in entry.info:
        return 0
    x = entry.info["ExcHet"]
    if isinstance(x, tuple):
        x = x[0]
    return x

def get_missing(entry):
    if "F_MISSING" not in entry.info:
        return 0
    x = entry.info["F_MISSING"]
    if isinstance(x, tuple):
        x = x[0]
    return x


def gt_filt(entry):
    """
    Returns True if Missingness is above .1 or number of samples is <= 1
    """
    n_present = 0
    for val in entry.samples.values():
        if 1 in val['GT']:
            n_present += 1
        if n_present > 1:
            return False
    return True

vcf = pysam.VariantFile(sys.argv[1])
out = pysam.VariantFile("/dev/stdout", 'w', header=vcf.header)

n_checked = 0
n_filt = 0
for entry in vcf:
    n_checked += 1
    if get_exchet(entry) < 0.05:
        n_filt += 1
        continue
    if get_missing(entry) >= .1:
        n_filt += 1
        continue
    sz = truvari.entry_size(entry)
    if sz < 50 or sz > 10_000:
        n_filt += 1
        continue
    if gt_filt(entry):
        n_filt += 1
        continue
    out.write(entry)
sys.stderr.write(f"check: {n_checked}; filt: {n_filt}\n")
