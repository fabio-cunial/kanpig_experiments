rm -rf bench_*
base=data/sv.vcf.gz
truvari bench -b $base --includebed data/HG002.f1_assembly_v2_genbank.dipcall.bed \
    --short --sizemin 50 --sizefilt 50 --sizemax 10000 -c data/HG002_resolved.vcf.gz -o bench_orig
truvari bench -b $base --includebed data/HG002.f1_assembly_v2_genbank.dipcall.bed \
    --short --sizemin 50 --sizefilt 50 --sizemax 10000 -c data/sniffles_regenotyped.vcf.gz -o bench_snif
truvari bench -b $base --includebed data/HG002.f1_assembly_v2_genbank.dipcall.bed \
    --short --sizemin 50 --sizefilt 50 --sizemax 10000 -c data/kanpig_regenotyped.vcf.gz -o bench_kanpig
truvari bench -b $base --includebed data/HG002.f1_assembly_v2_genbank.dipcall.bed \
    --short --sizemin 50 --sizefilt 50 --sizemax 10000 -c data/cutesv_regenotyped.vcf.gz -o bench_cutesv
truvari bench -b $base --includebed data/HG002.f1_assembly_v2_genbank.dipcall.bed \
    --short --sizemin 50 --sizefilt 50 --sizemax 10000 -c data/svjedi_regenotyped.vcf.gz -o bench_svjedi
