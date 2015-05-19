LOG_DIR=${HOME}"/projects/sequence_clustering/results/analysis_k562ctcf/logs/"
if [ ! -d "${LOG_DIR}" ]; then
    mkdir ${LOG_DIR}
fi

PROCESS_NAME=`echo "get_data"`
echo "/data/tools/R-3.1.1/bin/Rscript get_data.R" | \
qsub -l h_vmem=2g -v PATH -cwd -N ${PROCESS_NAME} \
        -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
		
