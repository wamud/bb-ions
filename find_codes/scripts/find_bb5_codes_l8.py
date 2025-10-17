import stim
import scipy
import pandas as pd
import itertools
import random
import sys
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



osd_options={ 'xyz_error_bias': [1, 1, 1], 'bp_method': "minimum_sum", 'ms_scaling_factor': 0.05, 'osd_method': "osd_cs", 'osd_order': 4, 'channel_update': None, 'seed': 42, 'max_iter': 9, 'tqdm_disable' : 1, 'error_bar_precision_cutoff': 1e-6 }


min_d_max = 8
min_k = 6

# for l in range(9, 11):
#     for m in range(5, 8):

l = 8
for m in range(5, 8):

    print(f"l = {l}, m = {m}")

    viable_entries = []
    count = 0

    temp_file = f"../found_codes/bb5_l{l}m{m}_partial.jsonl"
    results_file = f"../found_codes/bb5_l{l}m{m}.jsonl"

    # Supprimer ancien fichier temporaire si présent
    if os.path.exists(temp_file):
        os.remove(temp_file)

    all_ijs = list(itertools.product(range(l), range(m))) # contains every possible (i, j) for this l and m 
    random.seed(42)
    # random.shuffle(all_ijs)

    combos = itertools.product(all_ijs, all_ijs, all_ijs, all_ijs, all_ijs) 

    for combo in combos:

        ij0, ij1, ij2, ij3, ij4 = combo

        # Adding the same matrices together would change code stabiliser weight
        if ij0 == ij1:
            continue
        if ij2 == ij3 or ij2 == ij4 or ij3 == ij4:
            continue

        Aij = [ij0, ij1]
        Bij = [ij2, ij3, ij4]
        # (Could make Aij have three terms and Bij have two but if searching all the terms anyway that will just swap the roles of the left and right data qubits).

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

        if d_max < min_d_max:
            continue

        entry = {
            "nkd": [2 * l * m, k, bb5.min_logical_weight],
            "l" : l,
            "m" : m,
            "Aij": Aij,
            "Bij": Bij,
            # f"pL at p = {p}": bb5.osdw_logical_error_rate,
        }

        # print(entry)
        
        viable_entries.append(entry)
        count += 1

        if count % 10 == 0: # save every 10
            append_entries_to_json(viable_entries, temp_file)
            viable_entries = []
            print(f"{count} viable entries saved to {temp_file}")


    # Save any remaining entries after loop
    if viable_entries:
        append_entries_to_json(viable_entries, temp_file)

    # Rename temp file to final results file
    if os.path.exists(temp_file):
        os.rename(temp_file, results_file)
        print(f"Final results saved to {results_file}")
