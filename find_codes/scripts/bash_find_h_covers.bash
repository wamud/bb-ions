source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
conda activate base

for h in {2..10}; do
    python find_48_4_7_h_cover_codes.py "$h" &
done

for h in {2..10}; do
    python find_60_4_8_h_cover_codes.py "$h" &
done

wait