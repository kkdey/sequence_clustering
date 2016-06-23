
#create directory for logs (if it doesn't already exist)
LOG_DIR="logs/"
if [ ! -d "${LOG_DIR}" ]; then
	    mkdir ${LOG_DIR}
		fi

for i in `seq 1 100`; do
	    PROCESS_NAME=`echo "sim_clustering_$i"`
		    echo "/data/tools/R-3.1.1/bin/Rscript run_simple_simulation.R $i" | \
				qsub -l h_vmem=4g -v PATH -cwd -N ${PROCESS_NAME} \
				        -o ${LOG_DIR}${PROCESS_NAME}"_log.out" -e ${LOG_DIR}${PROCESS_NAME}"_log.err"
						done 

