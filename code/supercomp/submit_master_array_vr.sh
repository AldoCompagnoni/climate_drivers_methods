#!/bin/bash

# so test usage with an already pre-compiled stan model dir is:
# bash submit_master_array_vr.sh 1 1 /work/compagna/vr-4280462060

# to make a real run with its own compilation
# bash submit_master_array_vr.sh 1 35 airt log_lambda normal /work/compagna/vr-1519665825

NAME=vr

START_TASK=$1
END_TASK=$2
CLIMATE=$3
RESPONSE=$4
FAMILY=$5

if [[ -z $2 ]] ; then
  echo "usage: bash submit_master_array_vr.sh START_TASK END_TASK CLIM RESP FAM [MODEL_DIR]" >&2
  exit 1
fi

if [[ -n $6 ]] ; then
  MODEL_DIR="$6"
else
  # create unique output dir for this run
  MODEL_DIR="/work/$USER/$NAME-$(date +%s)"
  mkdir -p $MODEL_DIR

  # submit the stan compilation job and get its id
  compilation_job=$(qsub -terse submit_vr_compilation.sh $MODEL_DIR)
fi

# submit the array jobs which are waiting for the compilation job
qsub \
  ${compilation_job:+-hold_jid $compilation_job} \
  -t $START_TASK-$END_TASK \
  submit_array_vr_with_precompiled_models.sh "$MODEL_DIR" "$CLIMATE" "$RESPONSE" "$FAMILY"
