source /home/aforourk/Data/anaconda3/etc/profile.d/conda.sh
conda activate base

for h in {2..5}; do
    python find_h_cover_codes.py bb6_72_12_6 "$h" 16 18 &
done

for h in {6..10}; do
    python find_h_cover_codes.py bb6_72_12_6 "$h" 20 24 &
done


wait
