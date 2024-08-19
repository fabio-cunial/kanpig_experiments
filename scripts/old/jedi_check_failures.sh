#!/bin/bash
#
ADDRESS_FILE=$1
FAILURE_STRING="Killed"

#set -euxo #pipefail

cat ${ADDRESS_FILE} | tr '\t' ',' > tmp.txt
HEADER=""; IS_HEADER="1"
while read LINE; do
    if [ ${IS_HEADER} -eq 1 ]; then
        HEADER=${LINE}
        IS_HEADER="0"
    else
        SAMPLE=$(echo ${LINE} | cut -d , -f 1)
        for COLUMN in 2 3 4 5 6 7 8 9; do
            ADDRESS=$(echo ${LINE} | cut -d , -f ${COLUMN})
            if [ -n "${ADDRESS}" ]; then
                ADDRESS=$(dirname ${ADDRESS})
                ADDRESS=$(dirname ${ADDRESS})
                ADDRESS="${ADDRESS}/SvjediGenotyper.log"
                FAILED=$(gsutil cat ${ADDRESS} | grep ${FAILURE_STRING} | wc -l)
                if [ ${FAILED} -eq 1 ]; then
                    echo "Failure in sample=${SAMPLE} file=$(echo ${HEADER} | cut -d , -f ${COLUMN})"
                fi
            fi
        done
    fi
done < tmp.txt
