POP=$1
o_vcf=../filtered.1kgp.855.vcf.gz
k_vcf=../filtered.kanpig.855.vcf.gz
bed=../include.bed

bcftools view --force-samples -c 1 -S subset_${POP}.txt $k_vcf \
    | bcftools +fill-tags \
    | bcftools view -i "ExcHet >= 0.05 & MAF >= 0.01 & F_MISSING <= 0.1" -O z -o kanpig.${POP}.vcf.gz

tabix kanpig.${POP}.vcf.gz

python /users/u233287/fritz/english/code/kanpig_experiments/manuscript/1kgp/hwe_genotypes.py \
    kanpig.${POP}.vcf.gz \
    genotypes.kanpig.${POP}.jl  \
    $bed 

bcftools view --force-samples -c 1 -S subset_${POP}.txt $o_vcf \
    | bcftools +fill-tags \
    | bcftools view -i "ExcHet >= 0.05 & MAF >= 0.01 & F_MISSING <= 0.1" -O z -o 1kgp.${POP}.vcf.gz

tabix 1kgp.${POP}.vcf.gz

python /users/u233287/fritz/english/code/kanpig_experiments/manuscript/1kgp/hwe_genotypes.py \
    1kgp.${POP}.vcf.gz \
    genotypes.1kgp.${POP}.jl  \
    $bed 


