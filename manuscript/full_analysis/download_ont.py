import os
import pandas as pd
from google.cloud import storage

def download_blob(file_name, destination_file_name):
    """Downloads a blob from the bucket."""
    if not isinstance(file_name, str):
        return
    parts = file_name[5:].split('/')
    bucket_name = parts[0]
    source_blob_name = "/".join(parts[1:])

    # Initialize a client
    storage_client = storage.Client()

    # Get the bucket and the blob
    bucket = storage_client.bucket(bucket_name)
    blob = bucket.blob(source_blob_name)

    # Download the blob to a local file
    blob.download_to_filename(destination_file_name)

# manually have to add 32x for the 04 16 and 8x
FILES = ['02_dipcall_bed',
         '02_dipcall_sv_tbi',
         '02_dipcall_sv_vcf',
         '04_32x_sniffles_ccs_resolved_tbi',
         '04_32x_sniffles_ccs_resolved_vcf',
         '04_32x_sniffles_r10_resolved_tbi',
         '04_32x_sniffles_r10_resolved_vcf',
         '04_32x_sniffles_r10_resolved_tbi',
         '04_32x_sniffles_r10_resolved_vcf',
         '05_r10_cutesv_tbi',
         '05_r10_cutesv_vcf',
         '05_r10_jedi_tbi',
         '05_r10_jedi_vcf',
         '05_r10_kanpig_tbi',
         '05_r10_kanpig_vcf',
         '05_r10_sniffles_tbi',
         '05_r10_sniffles_vcf',
         '08_r10_cutesv_tbi',
         '08_r10_cutesv_vcf',
         '08_r10_jedi_tbi',
         '08_r10_jedi_vcf',
         '08_r10_kanpig_tbi',
         '08_r10_kanpig_vcf',
         '08_r10_sniffles_tbi',
         '08_r10_sniffles_vcf',
         '10_r10_cutesv_tbi',
         '10_r10_cutesv_vcf',
         '10_r10_jedi_tbi',
         '10_r10_jedi_vcf',
         '10_r10_kanpig_tbi',
         '10_r10_kanpig_vcf',
         '10_r10_sniffles_tbi',
         '10_r10_sniffles_vcf',

         '04_32x_sniffles_r9_resolved_tbi',
         '04_32x_sniffles_r9_resolved_vcf',
         '04_32x_sniffles_r9_resolved_tbi',
         '04_32x_sniffles_r9_resolved_vcf',
         '05_r9_cutesv_tbi',
         '05_r9_cutesv_vcf',
         '05_r9_jedi_tbi',
         '05_r9_jedi_vcf',
         '05_r9_kanpig_tbi',
         '05_r9_kanpig_vcf',
         '05_r9_sniffles_tbi',
         '05_r9_sniffles_vcf',
         '08_r9_cutesv_tbi',
         '08_r9_cutesv_vcf',
         '08_r9_jedi_tbi',
         '08_r9_jedi_vcf',
         '08_r9_kanpig_tbi',
         '08_r9_kanpig_vcf',
         '08_r9_sniffles_tbi',
         '08_r9_sniffles_vcf',
         '10_r9_cutesv_tbi',
         '10_r9_cutesv_vcf',
         '10_r9_jedi_tbi',
         '10_r9_jedi_vcf',
         '10_r9_kanpig_tbi',
         '10_r9_kanpig_vcf',
         '10_r9_sniffles_tbi',
         '10_r9_sniffles_vcf',
         ]

tmap = {'02_dipcall_bed': 'dipcall.bed',
        '02_dipcall_sv_tbi': 'dipcall.vcf.gz.tbi',
        '02_dipcall_sv_vcf': 'dipcall.vcf.gz'}
emap = {'vcf': 'vcf.gz',
        'tbi': 'vcf.gz.tbi'}

def name_resolver(name):
    """
    """
    if name.startswith("02"):
        return tmap[name]
    name = name.replace("_resolved", "").replace("_32x", "")
    section, tech, program, ext = name.split('_')
    ext = emap[ext]
    return f'{section}_{tech}_{program}.{ext}'

if __name__ == '__main__':
    data = pd.read_csv("paths_ont_8.txt", sep='\t')
    SAMPKEY = 'entity:card_08x_id'
    for i in FILES:
        assert i in data.columns, "missing %s" % (i)
    for _, row in data.iterrows():
        sample = row[SAMPKEY]
        if not os.path.exists(sample):
            os.mkdir(sample)
        for col in FILES:
            iname = row[col]
            oname = name_resolver(col)
            download_blob(iname, os.path.join(sample, oname))

