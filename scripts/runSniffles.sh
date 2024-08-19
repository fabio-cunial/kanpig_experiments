#!/bin/bash
#
N_ITERATIONS=$1
SELECTED_LOGICAL_CORES=$2
N_THREADS=$3
INPUT_VCF_GZ=$4
INPUT_BAM=$5
REFERENCE_FA=$6
OUTPUT_FILE=$7


for i in $(seq 1 ${N_ITERATIONS}); do
    ${TIME_COMMAND} taskset --cpu-list ${SELECTED_LOGICAL_CORES} sniffles --threads ${N_THREADS} --reference ${REFERENCE_FA} --input ${INPUT_BAM} --genotype-vcf ${INPUT_VCF_GZ} --vcf ${OUTPUT_FILE}
    rm -f ${OUTPUT_FILE}
done

