#!/bin/bash 
#set -beEuo pipefail

jobs=$1

cat $jobs | egrep -v '^#' | while read array exome gene anno ; do
    echo '#' $array $exome $gene $anno
    array_exome_consistency_check.sh $array $exome 2>/dev/null 
done

exit 0
array_exome_consistency_check.sh rs143435072 1:11193631:C:T 2>/dev/null
