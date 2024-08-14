base_bed=../giab_comp/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
hprc_bed=HG002.f1_assembly_v2_genbank.dipcall.mainchr.bed
bedtools genomecov -bga -g ~/code/references/grch38/GRCh38_1kg_mainchrs.fa.fai \
    -i <(cat $base_bed $hprc_bed  | bedtools sort) > coverage.bed
awk '$4==2' coverage.bed > dipcall_covered.bed
