#!/bin/bash

for ((l=6; l<=15; l++)); do
    for ((m=3; m<=l; m++)); do  # m <= l pour éviter les doublons symétriques
        
        n=$((2 * l * m))

        if (( n < 60 )); then
            continue  # skip combinations with n < 60
        elif (( n <= 80 )); then
            min_k=5
            min_d=9
        elif (( n <= 100 )); then
            min_k=6
            min_d=9
        elif (( n <= 120 )); then
            min_k=8
            min_d=9
        elif (( n <= 150 )); then
            min_k=9
            min_d=9
        else
            min_k=10
            min_d=10
        fi

        python3 find_bb5_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
