#!/bin/bash

for ((l=1; l<=8; l++)); do
    for ((m=1; m<l; m++)); do  # m <= l pour éviter les doublons symétriques

        n=2*l*m

        if (( n > 65 )); then
            continue
        elif (( n < 50)); then
            continue
        elif [[ $l -eq 1 && $m -eq 1 ]]; then # continue if m = l = 1
            continue
        elif (( m > l )); then # feel like there should be a symmetry argument that searching both (l=a, m=b) and (l=b, m=a) should be equivalent ? 
            continue
        fi

        min_k=4
        min_d=5

        python3.11 find_bb6_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
