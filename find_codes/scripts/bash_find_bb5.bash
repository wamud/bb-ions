#!/bin/bash

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <min_l> <max_l>"
    exit 1
fi

min_l=$1
max_l=$2

for ((l=min_l; l<=max_l; l++)); do
    for ((m=1; m<=l; m++)); do

        n=$((2 * l * m))

        if (( n > 65 )); then
            continue
        fi

        if (( n < 49 )); then
            continue
        else
            min_k=4
            min_d=6
        fi

        # for iHPC:
        source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
        conda activate base
        python find_bb5_codes.py "$l" "$m" "$min_k" "$min_d" &
        
        # for macbook air:
        python3.11 find_bb5_codes.py "$l" "$m" "$min_k" "$min_d" &
    done
done

wait