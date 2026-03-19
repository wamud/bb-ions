#!/bin/bash

for ((l=1; l<=3; l++)); do
    for ((m=1; m<=2; m++)); do  # m <= l pour éviter les doublons symétriques


        if [[ $l -eq 1 && $m -eq 1 ]]; then # continue if m = l = 1
            continue

        fi

        min_k=1
        min_d=2

        python3 find_bb6_codes.py "$l" "$m" "$min_k" "$min_d" &

    done
done

wait
