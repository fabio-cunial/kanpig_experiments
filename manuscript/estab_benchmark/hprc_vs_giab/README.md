
truvari bench --pick ac --includebed giab_comp/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed -b giab_comp/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz -c hprc.hg002.vcf.gz -o trubench

truvari refine -f ~/code/references/grch38/GRCh38_1kg_mainchrs.fa -r trubench/candidate.refine.bed -u -U -R -t 4 trubench

bedtools genomecov -bga -g ~/code/references/grch38/GRCh38_1kg_mainchrs.fa.fai -i <(cat giab_comp/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed HG002.f1_assembly_v2_genbank.dipcall.mainchr.bed | bedtools sort) > coverage.bed
# awk '$4==2' coverage.bed > dipcall_covered.bed

truvari bench --pick ac --includebed dipcall_covered.bed -b giab_comp/GRCh38_HG002-T2TQ100-V1.0_stvar.vcf.gz -c hprc.hg002.vcf.gz -o trubench

truvari refine -f ~/code/references/grch38/GRCh38_1kg_mainchrs.fa -r trubench/candidate.refine.bed -u -U -R -t 4 trubench

