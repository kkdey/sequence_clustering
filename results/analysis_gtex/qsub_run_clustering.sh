#!/bin/sh

MEM=40g
K=$1
sample_size=$2
region=$3

LOG_DIR="$HOME/projects/sequence_clustering/results/analysis_gtex/logs/"
if [ ! -d "${LOG_DIR}" ]; then
	mkdir ${LOG_DIR}
fi


PROCESS_NAME=`echo "run_clustering_$1_$3"`
echo "/data/tools/R-3.1.1/bin/Rscript run_clustering.R $1 $2 $3" | \
	qsub -l h_vmem=${MEM} -v PATH -cwd -N ${PROCESS_NAME} \
	    -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
