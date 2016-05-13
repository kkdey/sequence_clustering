#!/bin/sh

MEM=60g

LOG_DIR="$HOME/projects/sequence_clustering/results/analysis_gtex/logs/"
if [ ! -d "${LOG_DIR}" ]; then
	mkdir ${LOG_DIR}
fi

PROCESS_NAME=`echo "run_multiseq_pairs"`
echo "/data/tools/R-3.1.1/bin/Rscript run_multiseq_pairs_10_73576054_73611082.R" | \
	qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
	    -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
