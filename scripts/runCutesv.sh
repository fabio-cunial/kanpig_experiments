#!/bin/bash
#
N_ITERATIONS=$1
SELECTED_LOGICAL_CORES=$2
N_THREADS=$3
INPUT_VCF_GZ=$4
INPUT_BAM=$5
REFERENCE_FA=$6
OUTPUT_FILE=$7


CUTESV_PARAMS_CCS="--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --merge_ins_threshold 500 --merge_del_threshold 500"
MIN_SUPPORT="1"
SVLEN_MIN="30"
SVLEN_MAX="-1"
for i in $(seq 1 ${N_ITERATIONS}); do
    ${TIME_COMMAND} taskset --cpu-list ${SELECTED_LOGICAL_CORES} cuteSV --threads ${N_THREADS} -Ivcf ${INPUT_VCF_GZ} ${CUTESV_PARAMS_CCS} --genotype --min_support ${MIN_SUPPORT} --min_size ${SVLEN_MIN} --max_size ${SVLEN_MAX} ${INPUT_BAM} ${REFERENCE_FA} ${OUTPUT_FILE} ./cutesv_tmp
done
