### This is like find_bb6_codes.py except it searches for BB6 codes of the form of original Bravyi BB paper (i.e. not multivariate, A has terms x,y,y B has terms y,x,x (to various powers))


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
    print("This script requires inputs l, m, min_k, min_d")


l = int(sys.argv[1])
m = int(sys.argv[2])
min_k = int(sys.argv[3])
min_d_max = int(sys.argv[4])
unique_integer = os.getpid()

print(f"l = {l}, m = {m}")

filename = f"bb6_bravyi_l{l}_m{m}_k{min_k}_d{min_d_max}_id{unique_integer}"

temp_file = f"../found_codes/bravyi_BB6/{filename}_partial.jsonl"
results_file = f"../found_codes/bravyi_BB6/{filename}.jsonl"


# A has terms x,y,y B has terms y,x,x (to various powers)). So only need three values of each.
repeats = 3

ivalues = range(l)
ivectors = list(itertools.product(ivalues, repeat=repeats))
# random.shuffle(ivectors)

jvalues = range(m)
jvectors = list(itertools.product(jvalues, repeat=repeats))
# random.shuffle(jvectors)


# --- caches ---
seen_structures = set()   # (Aij, Bij) normalisés
seen_codes = set()        # matrices canoniques


# if j0 != 0 or i1 != 0 or i2 !=0 or i3 != 0 or j4 != 0 or j5 != 0:
    # continue
# A has terms x,y,y B has terms y,x,x (to various powers))
j0 = 0
i1 = 0
i2 = 0
i3 = 0
j4 = 0
j5 = 0

counting = 0
for loop in range(len(ivectors)):
    for count, jvec in enumerate(jvectors):
        ivec = ivectors[(count + loop) % len(ivectors)]

        i0, i4, i5 = ivec
        j1, j2, j3 = jvec


        # --- éliminer duplications internes ---
        if (i0 == i1 and j0 ==j1) or (i0 == i2 and j0 == j2) or (i1 == i2 and j1 == j2):
            continue
        if (i3 == i4 and j3 == j4) or (i3 == i5 and j3 == j5) or (i4 == i5 and j4 == j5):
            continue


        # # --- ordre lexicographique --- # skip any not in lexicographical order (avoids equivalent codes that are simply a re-ordering of terms)
        # if not ((i0, j0) < (i1, j1) < (i2, j2)):
        #     continue 
        # if not ((i3, j3) < (i4, j4) < (i5, j5)):
        #     continue 

        ## BB6 codes of the form of original BB paper (not multivariate, A is x,y,y B is y,x,x)   
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
