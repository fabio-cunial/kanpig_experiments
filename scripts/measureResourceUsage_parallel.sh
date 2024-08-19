#!/bin/bash
#
N_ITERATIONS="1"
TIME_COMMAND="/usr/bin/time --verbose"

INPUT_VCF_GZ="../practical/HG002_resolved.vcf.gz"
INPUT_BAM="../practical/HG002.bam"
REFERENCE_FA="../GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUTPUT_FILE="./outputfile"


# ----------------------------- KANPIG 0.3.1 -----------------------------------
rm -f ./kanpig_parallel_*.log
for N_THREADS in $(seq 64 -4 4) 2 1; do
    SELECTED_LOGICAL_CORES="$(( 64 - ${N_THREADS} ))-63"
    LOG_FILE="./kanpig_parallel_${N_THREADS}.log"
    rm -f ${LOG_FILE} ${OUTPUT_FILE}
    ${TIME_COMMAND} ./runKanpig.sh ${N_ITERATIONS} ${SELECTED_LOGICAL_CORES} ${N_THREADS} ${INPUT_VCF_GZ} ${INPUT_BAM} ${REFERENCE_FA} ${OUTPUT_FILE} &>> ${LOG_FILE}
    rm -f ${OUTPUT_FILE}
done
