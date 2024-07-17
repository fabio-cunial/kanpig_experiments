# Section 09. Presumably the same for Sections 06/07

gunzip -c sniffles_regenotyped.vcf.gz| cut -f1-10 | bgzip > sniffles_fixed.vcf.gz

grep -v "X\|Y\|_" HG002.f1_assembly_v2_genbank.dipcall.bed > dipcall.bed

truvari anno numneigh -i dipcall_bcftools_merge.vcf.gz | bgzip > dipcall_bcftools_merge_neigh.vcf.gz

truvari bench -b dipcall_bcftools_merge_neigh.vcf.gz \
	      -c sniffles_fixed.vcf.gz \
	      --pick ac -p 0.90 -P 0.90 --sizemin 50 --sizefilt 50 --sizemax 10000 \
	      -o bench --includebed dipcall.bed \
	      --no-ref a --bSample HG002

# For Section 5 I think we just have to drop the no-ref?
But we'll also need to set expectations from Section 4


# Note about the bcftools merge truth VCF

Of the 403,031 autosomal SVs between 50bp and 10kbp in the bcftools merge truth VCF, 

Of the 429,056 autosomal SVs between 50bp and 10kbp in the bcftools merge truth VCF, 
280,193 (65.3%) are at least 95% similar in terms of size and sequence of another SV of the same type and no further
than 500bp away. 

2024-07-11 10:15:46,412 [INFO] Zipped 429056 variants Counter({'base': 429056})
2024-07-11 10:15:46,412 [INFO] 68389 chunks of 429056 variants Counter({'base': 422090, '__filtered': 6966})
2024-07-11 10:15:46,421 [INFO] Wrote 227988 Variants
2024-07-11 10:15:46,421 [INFO] 201068 variants collapsed into 79125 variants
2024-07-11 10:15:46,421 [INFO] 396806 samples' FORMAT fields consolidated
2024-07-11 10:15:46,421 [INFO] Finished collapse
