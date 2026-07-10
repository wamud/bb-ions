source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
conda activate base

for h in {2..5}; do
    python find_h_cover_codes.py bb5_30_4_5 "$h" 10 10 &
done
for h in {6..10}; do
    python find_h_cover_codes.py bb6_72_12_6 "$h" 10 10 &
done

wait
