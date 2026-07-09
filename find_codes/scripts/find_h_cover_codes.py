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
from itertools import product

def factor_pairs(n): # finds pairs of factors of a number. E.g. factor_pairs(24) gives [(2, 12), (3, 8), (4, 6)]
    pairs = []
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            pairs.append((i, n // i))
    return pairs

def append_entries_to_json(entries, filename):
    """Append a list of dicts to a JSON lines file."""
    with open(filename, "a") as f:
        for entry in entries:
            # convert tuples to lists because JSON does not support tuples
            entry_clean = {k: list(v) if isinstance(v, tuple) else v for k, v in entry.items()}
            json.dump(entry_clean, f)
            f.write("\n")


def make_all_choices(terms, l, m, u, t):
    all_choices = []
    
    for i,j in terms:
        choices = []
    
        for r in range(u):
            i_prime = i + r * l
    
            for s in range(t):
                j_prime = j + s * m
                choices.append((i_prime, j_prime))
    
        all_choices.append(choices)
    
    return all_choices # gives a list of all possible choices for each term. E.g. all_choices[0] will give a list of all possible choices for the 0-th term


def all_lifts(Aij, Bij, og_l, og_m, u, t):
    A_choices = make_all_choices(Aij, og_l, og_m, u, t)
    B_choices = make_all_choices(Bij, og_l, og_m, u, t)
    
    for A_new in product(*A_choices):
        for B_new in product(*B_choices):
            yield list(A_new), list(B_new)

osd_options={ 'xyz_error_bias': [1, 1, 1], 'bp_method': "minimum_sum", 'ms_scaling_factor': 0.05, 'osd_method': "osd_cs", 'osd_order': 4, 'channel_update': None, 'seed': 42, 'max_iter': 9, 'tqdm_disable' : 1, 'error_bar_precision_cutoff': 1e-6}





code_name = str(sys.argv[1])
h = int(sys.argv[2])
min_k = int(sys.argv[3])
min_d = int(sys.argv[4])



if len(sys.argv) < 3:
    print("This script requires inputs code, h-cover, min_k, min_d. E.g. find_h_cover_codes.py gross_code 4 12 12 (see src/bb_ions/bbparamfuncs for possible code names)")
code = eval(f"{code_name}()")


print(f"Finding {h}-cover codes of [[{code.n}, {code.k}, {code.d_max}]] code with k ≥ {min_k}, d ≥ {min_d}")


og_l = code.l
og_m = code.m
og_Aij = code.Aij
og_Bij = code.Bij





filename = f"[[{code.n}, {code.k}, {code.d_max}]] {h}-covers min_k={min_k}, min_d={min_d}"
temp_file = f"../found_h_cover_codes/{filename}_partial.jsonl"
results_file = f"../found_h_cover_codes/{filename}.jsonl"
progress_file = f"../found_h_cover_codes/{filename}_progress.txt"

pairs = factor_pairs(h)
reversed_pairs = [(t, u) for u, t in pairs]
all_pairs = pairs + reversed_pairs
# print(all_pairs)

for u,t in all_pairs:

    l = u * og_l
    m = t * og_m

    print(f"l={l}, m={m}")


    seen_structures = set()   # (Aij, Bij) normalisés

    for Aij, Bij in all_lifts(og_Aij, og_Bij, og_l, og_m, u, t):
#        print(Aij)
#        print(Bij)


        # --- normalisation translation ---
        imin, jmin = Aij[0]

        Aij_norm = [((i - imin) % l, (j - jmin) % m) for (i, j) in Aij]
        Bij_norm = [((i - imin) % l, (j - jmin) % m) for (i, j) in Bij]

        Aij_norm = sorted(Aij_norm)
        Bij_norm = sorted(Bij_norm)
    
        # --- skip structures déjà vues ---
        key_struct = (tuple(Aij_norm), tuple(Bij_norm))
        if key_struct in seen_structures:
            continue
        seen_structures.add(key_struct)


        Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)

        k = 0

        if is_valid(Hx, Hz):
            k = find_k(Hx, Hz)

        if k < min_k:
            continue
        
        p = 0.1
        target_runs = 500

        with suppress_stdout():
            bb5 = css_decode_sim.css_decode_sim(hx = Hx, hz = Hz, error_rate = p, target_runs = target_runs, **osd_options)
            d_max = bb5.min_logical_weight

        if d_max < min_d:
            del bb5 
            gc.collect()
            continue

        entry = {
            "nkd": [2 * l * m, k, bb5.min_logical_weight],
            "l" : l,
            "m" : m,
            "Aij": Aij,
            "Bij": Bij,
            # f"pL at p = {p}": bb5.osdw_logical_error_rate,
        }

        print("  ", entry)

        append_entries_to_json([entry], temp_file)
        del bb5
        gc.collect()


# Rename temp file to final results file
if os.path.exists(temp_file):
    os.rename(temp_file, results_file)
    print(f"Final results saved to {results_file}")

