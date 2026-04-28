#!/bin/bash

for ((l=3; l<=3; l++)); do
    for ((m=3; m<=l; m++)); do  # m <= l pour éviter les doublons symétriques

        n=2*l*m



        if (( n != 18 )); then
            continue
        fi

        min_k=4
        min_d=2
        
        # source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
        #conda activate base
        # source ~/bb_env/bin/activate # for use on loan laptop

        python3 find_bb8_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
