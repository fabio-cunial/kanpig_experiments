#!/bin/bash
truvari bench -b sniffles.merge.callonly_fix.vcf.gz \
              -c sniffles.truvari.kanpig.callonly.vcf.gz \
              -o snif_kanpig --sizemin 50 --sizefilt 50 \
              --sizemax 10000 --includebed ../include.bed \
              --pctsize 0.9 --pctseq 0.9 --short
#!/bin/bash
truvari bench -b sniffles.merge.callonly_fix.vcf.gz \
              -c final-vcf.phased.vcf.gz \
              -o snif_1kgp --sizemin 50 --sizefilt 50 \
              --sizemax 10000 --includebed ../include.bed \
              --pctsize 0.9 --pctseq 0.9 --short
#!/bin/bash
truvari bench -b sniffles.truvari.kanpig.callonly.vcf.gz \
              -c final-vcf.phased.vcf.gz \
              -o kanpig_1kgp --sizemin 50 --sizefilt 50 \
              --sizemax 10000 --includebed ../include.bed \
              --pctsize 0.9 --pctseq 0.9 --short

