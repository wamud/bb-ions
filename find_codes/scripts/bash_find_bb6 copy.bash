#!/bin/bash

for ((l=21; l<=6; l++)); do
    for ((m=1; m<=6; m++)); do  # m <= l pour éviter les doublons symétriques


        if [[ $l -eq 1 && $m -eq 1 ]]; then # continue if m = l = 1
            continue
        elif (( 2*l*m > 48 )); then # Atm just searching for the BB6 codes used in the Delfosse BB5 paper
            continue
        elif (( m > l )); then # feel like there should be a symmetry argument that searching both (l=a, m=b) and (l=b, m=a) should be equivalent ? 
            continue
        fi

        min_k=4
        min_d=4

        python3.11 find_bb6_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
