#!/bin/bash
#
N_ITERATIONS="1"
SELECTED_LOGICAL_CORES="8-23"
N_THREADS="16"
TIME_COMMAND="/usr/bin/time --verbose"

INPUT_VCF_GZ="HG002_resolved.vcf.gz"
INPUT_BAM="HG002.bam"
INPUT_FASTQ="HG002.fastq"
REFERENCE_FA="../GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUTPUT_FILE="./outputfile"


# ---------------------------- SNIFFLES 2.3.3 ----------------------------------
LOG_FILE="./sniffles.log"

rm -f ${LOG_FILE} ${OUTPUT_FILE}
eval "$(conda shell.bash hook)"
conda activate sniffles
${TIME_COMMAND} ./runSniffles.sh ${N_ITERATIONS} ${SELECTED_LOGICAL_CORES} ${N_THREADS} ${INPUT_VCF_GZ} ${INPUT_BAM} ${REFERENCE_FA} ${OUTPUT_FILE} &>> ${LOG_FILE}
conda deactivate


# ------------------------- CUTESV (latest commit) -----------------------------
LOG_FILE="./cutesv.log"

mkdir ./cutesv_tmp
rm -f ${LOG_FILE} ${OUTPUT_FILE}
source ~/.kanpig/bin/activate
${TIME_COMMAND} ./runCutesv.sh ${N_ITERATIONS} ${SELECTED_LOGICAL_CORES} ${N_THREADS} ${INPUT_VCF_GZ} ${INPUT_BAM} ${REFERENCE_FA} ${OUTPUT_FILE} &>> ${LOG_FILE}
deactivate
rm -rf ${OUTPUT_FILE} ./cutesv_tmp


# ----------------------------- KANPIG 0.3.1 -----------------------------------
LOG_FILE="./kanpig.log"

rm -f ${LOG_FILE} ${OUTPUT_FILE}
${TIME_COMMAND} ./runKanpig.sh ${N_ITERATIONS} ${SELECTED_LOGICAL_CORES} ${N_THREADS} ${INPUT_VCF_GZ} ${INPUT_BAM} ${REFERENCE_FA} ${OUTPUT_FILE} &>> ${LOG_FILE}
rm -f ${OUTPUT_FILE}


# ------------------------------ JEDI 1.2.1 ------------------------------------
LOG_FILE="./jedi.log"

rm -f ${LOG_FILE} ${OUTPUT_FILE}
gunzip -c ${INPUT_VCF_GZ} > input.vcf
${TIME_COMMAND} ./runJedi.sh ${N_ITERATIONS} ${SELECTED_LOGICAL_CORES} ${N_THREADS} input.vcf ${INPUT_FASTQ} ${REFERENCE_FA} ${OUTPUT_FILE} &>> ${LOG_FILE}
rm -f ${OUTPUT_FILE} input.vcf
