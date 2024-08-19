#!/bin/bash
#


function printResources() {
    local LOG=$1
    local LABEL=$2
    
    TIME=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" ${LOG})
    TIME=${TIME:45:1000}
    TIME_S=$(echo ${TIME} | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }')
    RSS=$(grep "Maximum resident set size (kbytes):" ${LOG})
    RSS=${RSS:37:1000}
    
    echo "${LABEL},${TIME_S},${RSS}"
}


# Main program
for N_THREADS in $(seq 64 -4 4) 2 1; do
    printResources kanpig_parallel_${N_THREADS}.log ${N_THREADS}
done
