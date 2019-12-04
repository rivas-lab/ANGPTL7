#!/bin/bash
set -beEuo pipefail

SRCNAME=$(readlink -f $0)
PROGNAME=$(basename $SRCNAME)
VERSION="0.0.1"
NUM_POS_ARGS="0"

source "$(dirname ${SRCNAME})/GREML_misc.sh"

pheno=/oak/stanford/groups/mrivas/ukbb24983/phenotypedata/extras/iop/phe/INI2005254.phe
out=test
threads=4

gcta_reml $pheno $out $threads
