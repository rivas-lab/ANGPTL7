#!/bin/bash
#SBATCH    --job-name=GRM.combine
#SBATCH --output=logs_GRM/GRM.combine.%A.out
#SBATCH  --error=logs_GRM/GRM.combine.%A.err
#SBATCH  --nodes=1
#SBATCH  --cores=1
#SBATCH    --mem=8000
#SBATCH   --time=1:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail

############################################################
# config params
############################################################

pop="white_british"
n_parts=100
var_list="rs200058074_rs28991002_rs28991009_rs143435072"
# source "$(dirname ${SRCNAME})/GCTA_misc.sh"
source "/oak/stanford/groups/mrivas/users/ytanigaw/repos/rivas-lab/ANGPTL7/notebook/ukbb_GCTA/GCTA_misc.sh"

############################################################
# job start
############################################################

_SLURM_JOBID=${SLURM_JOBID:=0} # use 0 for default value (for debugging purpose)
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-start] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2

############################################################
# body
############################################################

part=$PI_SCRATCH/projects/ANGPTL7/ukbb_GCTA/grm_part/ukb24983_cal_cALL_v2_hg19.${pop}.${var_list}

gcta_grm_combine ${part} ${n_parts} ${GRM_WB_ANGPTL7}

############################################################
# job finish footer (for use with array-job module)
############################################################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [array-end] hostname=$(hostname) ; SLURM_JOBID=${_SLURM_JOBID}" >&2
