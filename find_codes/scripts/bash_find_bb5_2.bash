#!/bin/bash

# Nouvelles combinaisons restantes avec l > m
# (en respectant 50 < 2*l*m < 150)
declare -a pairs=(
    "13 2" "14 2" "15 2" "16 2" "17 2" "18 2" "19 2" "20 2" "21 2" "22 2" "23 2" "24 2" "25 2" "26 2" "27 2" "28 2" "29 2" "30 2" "31 2" "32 2" "33 2" "34 2" "35 2" "36 2" "37 2"
    "9 3"  "10 3" "11 3" "12 3" "13 3" "14 3" "15 3" "16 3" "17 3" "18 3" "19 3" "20 3" "21 3" "22 3" "23 3" "24 3" "25 3"
    "7 4"  "8 4"  "9 4"  "10 4" "11 4" "12 4" "13 4" "14 4" "15 4" "16 4" "17 4" "18 4"
    "6 5"  "7 5"  "8 5"  "9 5"  "10 5" "11 5"
)

for pair in "${pairs[@]}"; do
    read -r l m <<< "$pair"
    n=$((2 * l * m))

    if (( n < 60 )); then
        min_k=4
        min_d=8
    elif (( n >= 60 && n <= 80 )); then
        min_k=5
        min_d=9
    elif (( n > 80 && n <= 100 )); then
        min_k=6
        min_d=9
    elif (( n > 100 && n <= 120 )); then
        min_k=8
        min_d=9
    elif (( n > 120 && n <= 150 )); then
        min_k=9
        min_d=9
    else
        min_k=10
        min_d=10
    fi

    python3 find_bb5_codes.py "$l" "$m" "$min_k" "$min_d" &
done

wait

