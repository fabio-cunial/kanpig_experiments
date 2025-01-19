#!/bin/bash
#
set -x
WOMTOOL_PATH="/Users/fcunial/apps/cromwell/womtool-84.jar"

java -jar ${WOMTOOL_PATH} validate -l CutesvSingleSample.wdl
java -jar ${WOMTOOL_PATH} validate -l KanpigGenotyper.wdl
java -jar ${WOMTOOL_PATH} validate -l LrcallerGenotyper.wdl
java -jar ${WOMTOOL_PATH} validate -l SvjediGenotyper.wdl
java -jar ${WOMTOOL_PATH} validate -l CutesvGenotyper.wdl
java -jar ${WOMTOOL_PATH} validate -l SnifflesGenotyper.wdl
java -jar ${WOMTOOL_PATH} validate -l TruvariIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l SnifflesIntersample.wdl
java -jar ${WOMTOOL_PATH} validate -l Resolve.wdl
java -jar ${WOMTOOL_PATH} validate -l SnifflesSingleSample.wdl
