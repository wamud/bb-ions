#!/bin/bash

for ((l=2; l<=16; l++)); do
    for ((m=2; m<=4; m++)); do  # m <= l pour éviter les doublons symétriques

        n=2*l*m



        if (( n > 65 )); then
            continue
        elif (( n < 48)); then
            continue
        # elif (( m > l )); then # feel like there should be a symmetry argument that searching both (l=a, m=b) and (l=b, m=a) should be equivalent ? 
            # continue
        fi

        min_k=6
        min_d=8
	source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
	conda activate base
        python3 find_bb7_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
