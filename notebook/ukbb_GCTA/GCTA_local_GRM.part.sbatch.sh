#!/bin/bash
#SBATCH --job-name=GRM
#SBATCH --output=logs_GRM/GRM.%A_%a.out
#SBATCH  --error=logs_GRM/GRM.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=8G
#SBATCH --time=0:15:00
#SBATCH -p mrivas,normal,owners

set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="0"

############################################################
# config params
############################################################

ml load gcta
# source "$(dirname ${SRCNAME})/GCTA_misc.sh"
source "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ANGPTL7/notebook/ukbb_GCTA/GCTA_misc.sh"
pop="white_british"
var_list="rs200058074_rs28991002_rs28991009_rs143435072"
n_parts=100

cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
_SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID:=1}

############################################################
# job start
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname); SLURM_JOBID=${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${_SLURM_ARRAY_TASK_ID}" >&2

############################################################
# body
############################################################

grm=$PI_SCRATCH/projects/ANGPTL7/ukbb_GCTA/grm_part/ukb24983_cal_cALL_v2_hg19.${pop}.${var_list}
batch_num=${_SLURM_ARRAY_TASK_ID}

gcta_grm_part_extract ${GRM_CAL_BFILE} ${var_list} ${grm} ${pop} ${cores} ${batch_num} ${n_parts}

############################################################
# job finish footer (for use with array-job module)
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname); SLURM_JOBID=${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${_SLURM_ARRAY_TASK_ID}" >&2
