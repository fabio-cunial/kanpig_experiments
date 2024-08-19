#!/bin/bash
#
N_ITERATIONS=$1
SELECTED_LOGICAL_CORES=$2
N_THREADS=$3
INPUT_VCF_GZ=$4
INPUT_BAM=$5
REFERENCE_FA=$6
OUTPUT_FILE=$7


KANPIG_PARAMS_SINGLESAMPLE="--sizemin 20 --sizemax 10000 --chunksize 1000 --gpenalty 0.02 --hapsim 0.9999 --sizesim 0.90 --seqsim 0.85 --maxpaths 10000"
PLOIDY_BED="../grch38_male.bed"
for i in $(seq 1 ${N_ITERATIONS}); do
	taskset --cpu-list ${SELECTED_LOGICAL_CORES} ~/kanpig --threads ${N_THREADS} --ploidy-bed ${PLOIDY_BED} ${KANPIG_PARAMS_SINGLESAMPLE} --reference ${REFERENCE_FA} --input ${INPUT_VCF_GZ} --bam ${INPUT_BAM} --out ${OUTPUT_FILE}
done