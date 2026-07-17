#!/bin/bash

for ((m=2; m<=210; m++)); do
    for ((l=2; l<m+1; l++)); do  # m <= l pour éviter les doublons symétriques

        n=2*l*m



        if (( n >= 400 )); then
            continue
        elif (( n < 378)); then
            continue
        fi

        min_k=16
        min_d=18

        # On iHPC:
        source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
        conda activate base
        python find_bravyi_bb6_codes.py "$l" "$m" "$min_k" "$min_d"  &

        # # On laptop:
        # python3.11 find_bravyi_bb6_codes.py "$l" "$m" "$min_k" "$min_d"  &
    
    done
done


wait
