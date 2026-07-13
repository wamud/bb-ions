#!/bin/bash
for ((unique_integer=1; unique_integer<5; unique_integer++)); do
    for ((m=2; m<=210; m++)); do
        for ((l=2; l<m+1; l++)); do  # m <= l pour éviter les doublons symétriques

            n=2*l*m



            if (( n > 420 )); then
                continue
            elif (( n < 414)); then
                continue
            elif (( l > m )); then # feel like there should be a symmetry argument that searching both (l=a, m=b) and (l=b, m=a) should be equivalent ? 
                continue
            fi

            min_k=16
            min_d=18

            source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
            conda activate base

            python find_bb6_codes_2.py "$l" "$m" "$min_k" "$min_d" "$unique_integer"  &

        done
    done
done

wait
