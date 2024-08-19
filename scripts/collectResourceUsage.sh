#!/bin/bash
#
N_RUNS="5"
PBMM2_LOG="pbmm2.log"
SNIFFLES_LOG="sniffles.log"
CUTESV_LOG="cutesv.log"
JEDI_LOG="jedi.log"
KANPIG_LOG="kanpig.log"


function printResources() {
    local LOG=$1
    local LABEL=$2
    
    TIME=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" ${LOG})
    TIME=${TIME:45:1000}
    TIME_S=$(echo ${TIME} | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }')
    TIME_AVG=$( bc <<< "scale=2; ${TIME_S} / ${N_RUNS}" )
    RSS=$(grep "Maximum resident set size (kbytes):" ${LOG})
    RSS=${RSS:37:1000}
    
    echo "${LABEL},${TIME_AVG},${RSS}"
}


function printResources_minigraph() {
    local LOG=$1
    
    grep "Real time:" ${LOG} > minigraph.txt
    N_LINES=$(wc -l < minigraph.txt)
    MINIGRAPH_TIME_AVG="0"; MINIGRAPH_RAM_AVG="0"
    while read LINE; do
        LINE=${LINE:21:1000}
        MINIGRAPH_S=$(echo ${LINE} | cut -d ';' -f 1)
        MINIGRAPH_S=${MINIGRAPH_S:0:9}
        MINIGRAPH_RAM=$(echo ${LINE} | cut -d ';' -f 3)
        MINIGRAPH_RAM=${MINIGRAPH_RAM:11:6}
        MINIGRAPH_TIME_AVG=$( bc <<< "scale=2; ${MINIGRAPH_TIME_AVG} + ${MINIGRAPH_S}" )
        MINIGRAPH_RAM_AVG=$( bc <<< "scale=2; ${MINIGRAPH_RAM_AVG} + ${MINIGRAPH_RAM}" )
    done < minigraph.txt
    rm -f minigraph.txt
    MINIGRAPH_TIME_AVG=$( bc <<< "scale=2; ${MINIGRAPH_TIME_AVG} / ${N_LINES}" )
    MINIGRAPH_RAM_AVG=$( bc <<< "scale=2; (${MINIGRAPH_RAM_AVG} * 1000000) / ${N_LINES}" )

    echo "minigraph,${MINIGRAPH_TIME_AVG},${MINIGRAPH_RAM_AVG}"
}




# Main program
printResources ${PBMM2_LOG} pbmm2
printResources ${SNIFFLES_LOG} sniffles
printResources ${CUTESV_LOG} cutesv
printResources ${JEDI_LOG} jedi
printResources_minigraph ${JEDI_LOG}
printResources ${KANPIG_LOG} kanpig
