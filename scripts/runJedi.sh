#!/bin/bash
#
N_ITERATIONS=$1
SELECTED_LOGICAL_CORES=$2
N_THREADS=$3
INPUT_VCF=$4
INPUT_FASTQ=$5
REFERENCE_FA=$6
OUTPUT_FILE=$7

MIN_SUPPORT="1"
for i in $(seq 1 ${N_ITERATIONS}); do
    ${TIME_COMMAND} taskset --cpu-list ${SELECTED_LOGICAL_CORES} python ~/svjedigraph_patched/svjedi-graph.py --threads ${N_THREADS} --minsupport ${MIN_SUPPORT} --vcf ${INPUT_VCF} --ref ${REFERENCE_FA} --reads ${INPUT_FASTQ} --prefix ${OUTPUT_FILE}
done

