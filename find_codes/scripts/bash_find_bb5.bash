#!/bin/bash

for ((l=1; l<=3; l++)); do
    for ((m=1; m<=l; m++)); do  # m <= l pour éviter les doublons symétriques
        
        n=$((2 * l * m))

        if (( n > 12 )); then
            continue  
        elif (( n <= 12 )); then
            min_k=1
            min_d=2
        fi

        python3.11 find_bb5_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
