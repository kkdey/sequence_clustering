#!/bin/sh

MEM=20g
samplesize=$1
region=$2

LOG_DIR="$HOME/projects/sequence_clustering/data/gtex/logs/"
if [ ! -d "${LOG_DIR}" ]; then
	mkdir ${LOG_DIR}
fi

PROCESS_NAME=`echo "gtex_get_counts"`
echo "/data/tools/R-3.1.1/bin/Rscript get_counts_from_bam.R $samplesize $region" | \
	qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
    -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"






