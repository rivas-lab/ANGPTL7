#!/bin/bash
#SBATCH --job-name=HEreg
#SBATCH --output=logs/HEreg.local.%A_%a.out
#SBATCH  --error=logs/HEreg.local.%A_%a.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.2"
NUM_POS_ARGS="0"

############################################################
# config params
############################################################

ml load gcta
# source "$(dirname ${SRCNAME})/GCTA_misc.sh"
source "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ANGPTL7/notebook/ukbb_GCTA/GCTA_misc.sh"
pop="white_british"
var_list="rs200058074_rs28991002_rs28991009_rs143435072"

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

# pheno_file=/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/INI2005254.phe
# out=test.INI2005254

pheno_file="${GCTA_phe}"
pheno_idx="${_SLURM_ARRAY_TASK_ID}"
out="${DATA_D}/${pop}/$(get_phe_name ${pheno_idx}).${var_list}"

gcta_HEreg ${GRM_WB_ANGPTL7} $pheno_file $pheno_idx $out $cores

############################################################
# job finish footer (for use with array-job module)
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname); SLURM_JOBID=${_SLURM_JOBID}; SLURM_ARRAY_TASK_ID=${_SLURM_ARRAY_TASK_ID}" >&2
