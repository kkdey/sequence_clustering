#!/bin/sh

MEM=40g

LOG_DIR="$HOME/projects/sequence_clustering/results/analysis_gtex/logs/"
if [ ! -d "${LOG_DIR}" ]; then
	mkdir ${LOG_DIR}
fi


for i in ../../data/gtex/*.Robj; do
	PROCESS_NAME=`echo "run_multiseq_pairs"`
	echo "/data/tools/R-3.1.1/bin/Rscript run_multiseq_pairs.R $i" | \
		qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
		    -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
done
