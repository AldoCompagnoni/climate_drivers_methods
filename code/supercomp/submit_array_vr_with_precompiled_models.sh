#!/bin/bash 

#$  -N  array_vr
#$  -S  /bin/bash
#$  -l  h_rt=86400
#$  -l h_vmem=8G
#$  -cwd

# IMPORTAN: number of parallel environment needs to be SAME as binding linear number.
# first of the following two arguments specifies the number of cores to use.
#$  -binding linear:4
#$  -pe smp 4

#$ -o /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.log
#$ -e /work/$USER/$JOB_NAME-$JOB_ID-$TASK_ID.err

export LANG=en_US.UTF-8

if [[ -z $4 ]] ; then
  echo "usage: qsub submit_array_vr.sh /path/to/model_dir clim_var response family" >&2
  exit 1
fi

CODE_DIR="/home/$USER/sApropos_main"
OUTPUT_DIR="/work/$USER/$JOB_NAME-$JOB_ID"
# mkdir -p $OUTPUT_DIR (why doing this anyway?)

MODEL_DIR=$1
CLIMATE=$2
RESPONSE=$3
FAMILY=$4

echo "$@"

cd $MODEL_DIR

module load r/3.4.3-1

Rscript "$CODE_DIR/array_vr.R" \
	"$SGE_TASK_ID" \
	"/data/idiv_knight/sApropos/" \
	"$OUTPUT_DIR" \
	"$CLIMATE" \
	"$RESPONSE" \
        "$FAMILY" \
        "$CODE_DIR" \
	"${NSLOTS:-4}"

