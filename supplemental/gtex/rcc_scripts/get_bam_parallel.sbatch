#!/bin/sh

#SBATCH --time=36:00:00
#SBATCH --ntasks=32
#SBATCH --exclusive

module load parallel

# for large numbers of tasks the controlling node will have a large number
# of processes, so it will be necessary to change the user process limit
ulimit -u 10000

# the --exclusive to srun make srun use distinct CPUs for each job step
# -N1 -n1 allocates a single core to each task
srun="srun --exclusive -N1 -n1"

# --delay .2 prevents overloading the controlling node
# -j is the number of tasks parallel runs so we set it to $SLURM_NTASKS
# --joblog makes parallel create a log of tasks that it has already run
# --resume makes parallel use the joblog to resume from where it has left off
# the combination of --joblog and --resume allow jobs to be resubmitted if
# necessary and continue from where they left off
parallel="parallel --delay .2 -j $SLURM_NTASKS --joblog get_bam.log --resume"

# this runs the parallel command we want
# in this case, we are running a script named runtask
# parallel uses ::: to separate options. Here {1..128} is a shell expansion
# so parallel will run the runtask script for the numbers 1 through 128
# {1} is the first argument
# as an example, the first job will be run like this:
# srun --exclusive -N1 -n1 ./runtask arg1:1 > runtask.1
$parallel "for i in `ls jobs-*`; do sbatch $i; done"
