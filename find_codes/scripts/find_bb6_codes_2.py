# import stim
# import scipy
# import pandas as pd
# import itertools
# import random
# import sys
# import gc
# import os
# import json
# sys.path.append(os.path.abspath("../../src"))
# from bb_ions import *



# def append_entries_to_json(entries, filename):
#     """Append a list of dicts to a JSON lines file."""
#     with open(filename, "a") as f:
#         for entry in entries:
#             # convert tuples to lists because JSON does not support tuples
#             entry_clean = {k: list(v) if isinstance(v, tuple) else v for k, v in entry.items()}
#             json.dump(entry_clean, f)
#             f.write("\n")


# osd_options={ 
#     'xyz_error_bias': [1, 1, 1], 
#     'bp_method': "minimum_sum", 
#     'ms_scaling_factor': 0.05, 
#     'osd_method': "osd_cs", 
#     'osd_order': 4, 
#     'channel_update': None, 
#     'seed': 42, 
#     'max_iter': 9, 
#     'tqdm_disable' : 1, 
#     'error_bar_precision_cutoff': 1e-6
#     }


# if len(sys.argv) < 4:
#     print("This script requires inputs l, m, min_k, min_d")


# l = int(sys.argv[1])
# m = int(sys.argv[2])
# min_k = int(sys.argv[3])
# min_d_max = int(sys.argv[4])


# print(f"l = {l}, m = {m}")


# filename = f"bb6_l{l}_m{m}_k{min_k}_d{min_d_max}"

# temp_file = f"../found_codes/{filename}_partial.jsonl"
# results_file = f"../found_codes/{filename}.jsonl"
# progress_file = f"../found_codes/l{l}_m{m}_progress.txt"



# # We will make ivectors containing i0, i1, i2, i3, i4, i5 and jvectors containing j0, j1, j2, j3, j4, j5
# # which are the powers of the terms in the matrices A and B, i.e. 
# # A = x^i0 * y^j0 + x^i1 * y^j1 + x^i2 * y^j2
# # B = x^i3 * y^j3 + x^i4 * y^j4 + x^i5 * y^j5
# # In turn, A and B will be used in the Bicycle Bivariate code's parity check matrices as Hx = [A|B] and Hz = [B^T|A^T]

# random.seed(42)

# ivalues = range(l)
# ivectors = list(itertools.product(ivalues, repeat = 6))
# random.shuffle(ivectors)

# jvalues = range(m)
# jvectors = list(itertools.product(jvalues, repeat = 6))
# random.shuffle(jvectors)

# seen_structures = set()

# # To test all the possible combinations but to do so in a random order (without loading the full list of combos and using random.shuffle() as its size is prohibitive) lets do this funky nested for loop:
# for loop in range(len(ivectors)):
#     for count, jvec in enumerate(jvectors):
#         ivec = ivectors[(count + loop) % len(ivectors)]

#         if count % 100 == 0:
#             with open(progress_file, "w") as f:
#                 f.write(f"loop = {loop}, count = {count}")


#         i0 = ivec[0]
#         i1 = ivec[1]
#         i2 = ivec[2]
#         i3 = ivec[3]
#         i4 = ivec[4]
#         i5 = ivec[5]
        
#         j0 = jvec[0]
#         j1 = jvec[1]
#         j2 = jvec[2]
#         j3 = jvec[3]
#         j4 = jvec[4]
#         j5 = jvec[5]

#         # Skip values where the same term appears twice in the same matrix (this would change the weight of the stabilisers as when they're added together mod 2 their ones cancel out in the parity check matrices)
#         if (i0, j0) == (i1, j1) or (i0, j0) == (i2, j2) or (i1, j1) == (i2, j2):
#             continue
#         if (i3, j3) == (i4, j4) or (i3, j3) == (i5, j5) or (i4, j4) == (i5, j5):
#             continue

#         # Also equivalent polynomials with reordered terms need not be repeated. I.e. x^1y^2 + x^3 is equivalent to x^3 + x^1y^2.
#         # To avoid repeats let's just search terms that are in ascending order lexicographically (compare first element of tuple, if equal then compare second element of tuple)
#         if not ((i0, j0) < (i1, j1) < (i2, j2)):
#             continue 
#         if not ((i3, j3) < (i4, j4) < (i5, j5)):
#             continue 

#         Aij = [(i0, j0), (i1, j1), (i2, j2)]
#         Bij = [(i3, j3), (i4, j4), (i5, j5)]

#         # --- normalisation (translation globale) ---
#         imin, jmin = Aij[0]

#         Aij = [((i - imin) % l, (j - jmin) % m) for (i, j) in Aij]
#         Bij = [((i - imin) % l, (j - jmin) % m) for (i, j) in Bij]

#         Aij = sorted(Aij)
#         Bij = sorted(Bij)

#         key = (tuple(Aij), tuple(Bij))

#         if key in seen_structures:
#             continue

#         seen_structures.add(key)


#         Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)

#         k = 0

#         if is_valid(Hx, Hz):
#             k = find_k(Hx, Hz)

#         if k < min_k:
#             continue
        
#         p = 0.1
#         target_runs = 500

#         with suppress_stdout():

#             bb6 = css_decode_sim.css_decode_sim(hx = Hx, hz = Hz, error_rate = p, target_runs = target_runs, **osd_options)
#             d_max = bb6.min_logical_weight

#         if d_max < min_d_max:
#             del bb6 
#             gc.collect()
#             continue

#         n = 2 * l * m
#         if n == 18 and k == 4 and d_max == 2:
#             del bb6 
#             gc.collect()
#             continue
#         if n == 18 and k == 4 and d_max == 4:
#             del bb6 
#             gc.collect()
#             continue
#         if n == 18 and k == 12 and d_max == 2:
#             del bb6 
#             gc.collect()
#             continue
        
        

#         entry = {
#             "nkd": [2 * l * m, k, bb6.min_logical_weight],
#             "l" : l,
#             "m" : m,
#             "Aij": Aij,
#             "Bij": Bij,
#             # f"pL at p = {p}": bb6.osdw_logical_error_rate,
#         }

#         print(entry)

#         append_entries_to_json([entry], temp_file)
#         del bb6
#         gc.collect()


# # Rename temp file to final results file
# if os.path.exists(temp_file):
#     os.rename(temp_file, results_file)
#     print(f"Final results saved to {results_file}")



import stim
import scipy
import pandas as pd
import itertools
import random
import sys
import gc
import os
import json
import numpy as np

sys.path.append(os.path.abspath("../../src"))
from bb_ions import *


def append_entries_to_json(entries, filename):
    with open(filename, "a") as f:
        for entry in entries:
            entry_clean = {k: list(v) if isinstance(v, tuple) else v for k, v in entry.items()}
            json.dump(entry_clean, f)
            f.write("\n")


def canonical_matrix(H, max_iter=3):
    for _ in range(max_iter):
        # sort rows
        H = H[np.lexsort(H.T[::-1])]
        # sort columns
        H = H[:, np.lexsort(H[::-1])]
    return H


osd_options = { 
    'xyz_error_bias': [1, 1, 1], 
    'bp_method': "minimum_sum", 
    'ms_scaling_factor': 0.05, 
    'osd_method': "osd_cs", 
    'osd_order': 4, 
    'channel_update': None, 
    'seed': 42, 
    'max_iter': 100, 
    'tqdm_disable': 1, 
    'error_bar_precision_cutoff': 1e-6
}


if len(sys.argv) < 5:
    print("This script requires inputs l, m, min_k, min_d, unique_integer")


l = int(sys.argv[1])
m = int(sys.argv[2])
min_k = int(sys.argv[3])
min_d_max = int(sys.argv[4])
index = int(sys.argv[5])

print(f"l = {l}, m = {m}")

filename = f"bb6_{index}_l{l}_m{m}_k{min_k}_d{min_d_max}"

temp_file = f"../found_codes/around_416/{filename}_partial.jsonl"
results_file = f"../found_codes/around_416/{filename}.jsonl"


ivalues = range(l)
ivectors = list(itertools.product(ivalues, repeat=6))
random.shuffle(ivectors)

jvalues = range(m)
jvectors = list(itertools.product(jvalues, repeat=6))
random.shuffle(jvectors)


# --- caches ---
seen_structures = set()   # (Aij, Bij) normalisés
seen_codes = set()        # matrices canoniques


for loop in range(len(ivectors)):
    for count, jvec in enumerate(jvectors):
        ivec = ivectors[(count + loop) % len(ivectors)]



        i0, i1, i2, i3, i4, i5 = ivec
        j0, j1, j2, j3, j4, j5 = jvec


        # --- éliminer duplications internes ---
        if (i0, j0) == (i1, j1) or (i0, j0) == (i2, j2) or (i1, j1) == (i2, j2):
            continue
        if (i3, j3) == (i4, j4) or (i3, j3) == (i5, j5) or (i4, j4) == (i5, j5):
            continue


        # --- ordre lexicographique ---
        if not ((i0, j0) < (i1, j1) < (i2, j2)):
            continue 
        if not ((i3, j3) < (i4, j4) < (i5, j5)):
            continue 

        # BB6 codes of the form of original BB paper (not multivariate, A is x,y,y B is y,x,x)    
        # if j0 != 0 or i1 != 0 or i2 !=0 or i3 != 0 or j4 != 0 or j5 != 0:
            # continue


        Aij = [(i0, j0), (i1, j1), (i2, j2)]
        Bij = [(i3, j3), (i4, j4), (i5, j5)]


        # --- normalisation translation ---
        imin, jmin = Aij[0]

        Aij = [((i - imin) % l, (j - jmin) % m) for (i, j) in Aij]
        Bij = [((i - imin) % l, (j - jmin) % m) for (i, j) in Bij]

        Aij = sorted(Aij)
        Bij = sorted(Bij)


        # --- skip structures déjà vues ---
        key_struct = (tuple(Aij), tuple(Bij))
        if key_struct in seen_structures:
            continue
        seen_structures.add(key_struct)


        # --- construire matrices ---
        Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)


        # --- canonicalisation (lignes + colonnes) ---
        Hx_c = canonical_matrix(Hx)
        Hz_c = canonical_matrix(Hz)

        key_code = (Hx_c.tobytes(), Hz_c.tobytes())

        if key_code in seen_codes:
            continue
        seen_codes.add(key_code)


        # --- calcul k ---
        k = 0
        if is_valid(Hx, Hz):
            k = find_k(Hx, Hz)

        if k < min_k:
            continue


        # --- distance ---
        p = 0.1
        target_runs = 1000

        with suppress_stdout():
            bb6 = css_decode_sim.css_decode_sim(
                hx=Hx, hz=Hz,
                error_rate=p,
                target_runs=target_runs,
                **osd_options
            )
            d_max = bb6.min_logical_weight


        if d_max < min_d_max:
            del bb6
            gc.collect()
            continue


        n = 2 * l * m



        entry = {
            "nkd": [n, k, d_max],
            "l": l,
            "m": m,
            "Aij": Aij,
            "Bij": Bij,
        }

        print(entry)

        append_entries_to_json([entry], temp_file)

        del bb6
        gc.collect()


if os.path.exists(temp_file):
    os.rename(temp_file, results_file)
    print(f"Final results saved to {results_file}")
