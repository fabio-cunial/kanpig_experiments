#!/bin/bash
#
REMOTE_ROOT=$1  # E.g. "gs://.../submissions/.../SvjediGenotyper"
FAILURE_STRING="Killed"
END_STRING="delocalization"

set -euo pipefail


N_FAILED="0"; N_RUNNING="0"; N_NO_LOG="0"
gsutil ls ${REMOTE_ROOT} > jobs.txt
while read LINE; do
    #FOUND=$(gsutil cp ${LINE}call-SvjediGenotyper/SvjediGenotyper.log ./log.txt && echo 1 || echo 0)    
    FOUND=$(gsutil cp ${LINE}call-SvjediGenotyperImpl/SvjediGenotyperImpl.log ./log.txt && echo 1 || echo 0)    
    if [ ${FOUND} = "1" ]; then
        ENDED=$(grep ${END_STRING} log.txt && echo 1 || echo 0)
        if [ "${ENDED}" = "0" ]; then
            echo "          Still running"
            N_RUNNING=$(( ${N_RUNNING} + 1 ))
        else
            FAILED=$(grep ${FAILURE_STRING} log.txt && echo 1 || echo 0)
            if [ "${FAILED}" = "0" ]; then
                echo "Success"
            else
                echo "          FAILED"
                N_FAILED=$(( ${N_FAILED} + 1 ))
            fi
        fi
    else
        echo "          The job has no log"
        N_NO_LOG=$(( ${N_NO_LOG} + 1 ))
    fi
    rm -f log.txt
done < jobs.txt
rm -f jobs.txt
echo ""
echo "${N_FAILED} jobs have failed"
echo "${N_RUNNING} jobs are still running"
echo "${N_NO_LOG} jobs have no log"
