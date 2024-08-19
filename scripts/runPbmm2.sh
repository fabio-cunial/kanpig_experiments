#!/bin/bash
#
N_ITERATIONS=$1
SELECTED_LOGICAL_CORES=$2
N_THREADS=$3
INPUT_FASTQ=$4
REFERENCE_FA=$5
OUTPUT_FILE=$6


for i in $(seq 1 ${N_ITERATIONS}); do
    ${TIME_COMMAND} taskset --cpu-list ${SELECTED_LOGICAL_CORES} pbmm2 align --num-threads ${N_THREADS} --preset CCS --sort --sample sample ${REFERENCE_FA} ${INPUT_FASTQ} ${OUTPUT_FILE}
    #${TIME_COMMAND} taskset --cpu-list ${SELECTED_LOGICAL_CORES} samtools calmd -@ ${N_THREADS} --no-PG -b ${OUTPUT_FILE} ${REFERENCE_FA} > tmp.bam
    ${TIME_COMMAND} taskset --cpu-list ${SELECTED_LOGICAL_CORES} samtools index -@ ${N_THREADS} ${OUTPUT_FILE}
    rm -f ${OUTPUT_FILE} tmp.bam*
done

