#!/bin/bash

for ((l=6; l<=12; l++)); do
    for ((m=3; m<=12; m++)); do
        n=$((2 * l * m))

        if (( n < 60 )); then
            min_k=4
            min_d=7
        elif (( n >= 60 && n <= 80 )); then
            min_k=10
            min_d=6
        elif (( n > 80 && n <= 120 )); then
            min_k=8
            min_d=8
        elif (( n > 120 && n <= 150 )); then
            min_k=8
            min_d=10
        else
            min_k=10
            min_d=12
        fi

        python3 find_bb5_codes.py "$l" "$m" "$min_k" "$min_d" &
    done
done

wait
