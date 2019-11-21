#!/bin/bash
#SBATCH --job-name=IOP_GWAS
#SBATCH --output=log/IOP_GWAS.%A.out
#SBATCH  --error=log/IOP_GWAS.%A.err
#SBATCH --nodes=1
#SBATCH --cores=4
#SBATCH --mem=32000
#SBATCH --time=1-0:00:00
#SBATCH -p mrivas,normal

set -beEuo pipefail
cores=$( cat $0 | egrep '^#SBATCH --cores='  | awk -v FS='=' '{print $NF}' )
mem=$(   cat $0 | egrep '^#SBATCH --mem='    | awk -v FS='=' '{print $NF}' )

echo "[$0 $(date +%Y%m%d-%H%M%S)] [start] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}" >&2
#################
#gwas_py=/home/users/ytanigaw/repos/rivas-lab/ukbb-tools/04_gwas/gwas.py
#GBE_ID=INI5254

GBE_ID=$1
if [ $# -gt 1 ] ; then pop=$2 ; else pop=white_british ; fi

for runtype in array exome exome-gatk ; do
for pop in e_asian s_asian african non_british_white white_british ; do

if [ ! -d $(readlink -f $(pwd)/out/${pop}) ] ; then mkdir -p $(readlink -f $(pwd)/out/${pop}) ; fi

${gwas_py} \
    --run-${runtype} --run-now --memory $mem --cores $cores \
    --pheno $(readlink -f $(pwd)/phe/${GBE_ID}.phe) \
    --out $(readlink -f $(pwd)/out/${pop})  \
    --population ${pop} \
    --log-dir $(readlink -f $(pwd)/log)
done
done

#################
echo "[$0 $(date +%Y%m%d-%H%M%S)] [end] hostname = $(hostname) SLURM_JOBID = ${SLURM_JOBID:=0}" >&2

