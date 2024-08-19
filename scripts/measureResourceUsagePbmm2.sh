#!/bin/bash
#
N_ITERATIONS="5"
SELECTED_LOGICAL_CORES="4"
N_THREADS="1"
TIME_COMMAND="/usr/bin/time --verbose"

INPUT_FASTQ="HG002_chr1.fastq"
REFERENCE_FA="chr1.fa"
OUTPUT_FILE="./outputfile.bam"


# ---------------------------- PBMM2 1.13.0 ------------------------------------
LOG_FILE="./pbmm2.log"

rm -f ${LOG_FILE} ${OUTPUT_FILE}
eval "$(conda shell.bash hook)"
conda activate pbmm2
${TIME_COMMAND} ./runPbmm2.sh ${N_ITERATIONS} ${SELECTED_LOGICAL_CORES} ${N_THREADS} ${INPUT_FASTQ} ${REFERENCE_FA} ${OUTPUT_FILE} &>> ${LOG_FILE}
conda deactivate
