#!/bin/bash 

#$  -N  sim_regular
#$  -S  /bin/bash
#$  -l  h_rt=900
#$  -l h_vmem=32G
#$  -cwd

# IMPORTANT: number of parallel environment needs to be SAME as binding linear number.
# first of the following two arguments specifies the number of cores to use.
#$  -pe smp 4
#$  -binding linear:4

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -e /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.err

module load foss/2018b R/3.6.0

Rscript /home/compagna/sApropos_main/sim_regular/sim_regular_eve.R \
	"$SGE_TASK_ID" \
	"/gpfs1/work/$USER/sim_regular/$JOB_NAME-$JOB_ID"
