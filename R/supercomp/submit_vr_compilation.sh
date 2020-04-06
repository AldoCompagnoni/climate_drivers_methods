#!/bin/bash

#$  -N  vr_compilation
#$  -S  /bin/bash
#$  -l  h_rt=86400
#$  -l h_vmem=32G,highmem
#$  -cwd

# IMPORTANT: number of parallel environment needs to be SAME as binding linear number.
# first of the following two arguments specifies the number of cores to use.
#$  -binding linear:1

#$ -o /work/$USER/$JOB_NAME-$JOB_ID.log
#$ -e /work/$USER/$JOB_NAME-$JOB_ID.err

export LANG=en_US.UTF-8

if [[ -z $1 ]] ; then
  echo "usage: qsub submit_vr_compilation.sh /path/to/dir" >&2
  exit 1
fi

CODE_DIR=/home/$USER/sApropos_main

MODEL_DIR=$1

cp $CODE_DIR/*.stan $MODEL_DIR

cd $MODEL_DIR

module load r/3.4.3-1

Rscript $CODE_DIR/compilation_vr.R
