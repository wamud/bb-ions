import stim
import scipy
import pandas as pd
import itertools
import random
import sys
import gc
import os
import json
sys.path.append(os.path.abspath("../../src"))
from bb_ions import *



def append_entries_to_json(entries, filename):
    """Append a list of dicts to a JSON lines file."""
    with open(filename, "a") as f:
        for entry in entries:
            # convert tuples to lists because JSON does not support tuples
            entry_clean = {k: list(v) if isinstance(v, tuple) else v for k, v in entry.items()}
            json.dump(entry_clean, f)
            f.write("\n")


osd_options={ 
    'xyz_error_bias': [1, 1, 1], 
    'bp_method': "minimum_sum", 
    'ms_scaling_factor': 0.05, 
    'osd_method': "osd_cs", 
    'osd_order': 4, 
    'channel_update': None, 
    'seed': 42, 
    'max_iter': 9, 
    'tqdm_disable' : 1, 
    'error_bar_precision_cutoff': 1e-6
    }


if len(sys.argv) < 4:
    print("This script requires inputs l, m, min_k, min_d")


l = int(sys.argv[1])
m = int(sys.argv[2])
min_k = int(sys.argv[3])
min_d_max = int(sys.argv[4])


print(f"l = {l}, m = {m}")


filename = f"bb6_l{l}_m{m}_k{min_k}_d{min_d_max}"

temp_file = f"../found_codes/{filename}_partial.jsonl"
results_file = f"../found_codes/{filename}.jsonl"
progress_file = f"../found_codes/l{l}_m{m}_progress.txt"



# We will make ivectors containing i0, i1, i2, i3, i4, i5 and jvectors containing j0, j1, j2, j3, j4, j5
# which are the powers of the terms in the matrices A and B, i.e. 
# A = x^i0 * y^j0 + x^i1 * y^j1 + x^i2 * y^j2
# B = x^i3 * y^j3 + x^i4 * y^j4 + x^i5 * y^j5
# In turn, A and B will be used in the Bicycle Bivariate code's parity check matrices as Hx = [A|B] and Hz = [B^T|A^T]

# cover code values:
i0_og, i1_og, i2_og, i3_og, i4_og, i5_og = 0, 2, 2, 5, 9, 11  # the x-power values of the original code 
j0_og, j1_og, j2_og, j3_og, j4_og, j5_og = 0, 8, 6, 0, 9, 2

# --- caches ---
seen_structures = set()   # (Aij, Bij) normalisés
seen_codes = set()        # matrices canoniques


for b in [0, 1]:
    for c in [0, 1]:
        for d in [0, 1]:
            j0, j1, j2, j3, j4, j5 = j0_, 10 - b * 9 , 17 - c * 9 , 5 - d * 9 , 0, 0

            # Skip values where the same term appears twice in the same matrix (this would change the weight of the stabilisers as when they're added together mod 2 their ones cancel out in the parity check matrices)
            if (i0, j0) == (i1, j1) or (i0, j0) == (i2, j2) or (i1, j1) == (i2, j2):
                continue
            if (i3, j3) == (i4, j4) or (i3, j3) == (i5, j5) or (i4, j4) == (i5, j5):
                continue

            ## Commented the below lines out becasue original polynomials were not lexicographically ordered
            # # Equivalent polynomials with reordered terms need not be repeated. I.e. x^1y^2 + x^3 is equivalent to x^3 + x^1y^2.
            # # To avoid repeats let's just search terms that are in ascending order lexicographically (compare first element of tuple, if equal then compare second element of tuple)
            # if not ((i0, j0) < (i1, j1) < (i2, j2)):
            #     continue 
            # if not ((i3, j3) < (i4, j4) < (i5, j5)):
            #     continue 

            Aij_orig = [(i0, j0), (i1, j1), (i2, j2)] #pre-normalisation values
            Bij_orig = [(i3, j3), (i4, j4), (i5, j5)]


            # --- normalisation translation for avoidind duplicates ---
            # Translates all i and j values so that Aij[0] = (0, 0) -- this translation is equivalent to just reordering the columns. We do this because then we can compare codes. When we save the entry however we will save the original Aij Bij values rather than these renormalised ones which can look weird. (I.e. you might have identity in the first term but then very high powers in other terms )

            imin, jmin = Aij_orig[0]

            Aij_norm = [((i - imin) % l, (j - jmin) % m) for (i, j) in Aij_orig]
            Bij_norm = [((i - imin) % l, (j - jmin) % m) for (i, j) in Bij_orig]

            Aij_norm = sorted(Aij_norm)
            Bij_norm = sorted(Bij_norm)


            # --- skip structures déjà vues ---
            key_struct = (tuple(Aij_norm), tuple(Bij_norm))
            if key_struct in seen_structures:
                continue
            seen_structures.add(key_struct)


            print(f"Aij = {Aij_orig}")
            print(f"Bij = {Bij_orig}")


            Hx, Hz = make_parity_check_matrices(l, m, Aij_orig, Bij_orig)


            k = 0

            if is_valid(Hx, Hz):
                k = find_k(Hx, Hz)

            if k < min_k:
                continue
            
            p = 0.1
            target_runs = 500

            with suppress_stdout():

                bb6 = css_decode_sim.css_decode_sim(hx = Hx, hz = Hz, error_rate = p, target_runs = target_runs, **osd_options)
                d_max = bb6.min_logical_weight

            if d_max < min_d_max:
                del bb6 
                gc.collect()
                continue

            # n = 2 * l * m
            # if n == 18 and k == 4 and d_max == 2:
            #     del bb6 
            #     gc.collect()
            #     continue


            entry = {
                "nkd": [2 * l * m, k, bb6.min_logical_weight],
                "l" : l,
                "m" : m,
                "Aij": Aij_orig,
                "Bij": Bij_orig,
                # f"pL at p = {p}": bb6.osdw_logical_error_rate,
            }

            print(entry)

            append_entries_to_json([entry], temp_file)
            del bb6
            gc.collect()


# Rename temp file to final results file
if os.path.exists(temp_file):
    os.rename(temp_file, results_file)
    print(f"Final results saved to {results_file}")

