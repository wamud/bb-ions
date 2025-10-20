# import stim
# import scipy
# import pandas as pd
# import itertools
# import random
# import sys
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



# osd_options={ 'xyz_error_bias': [1, 1, 1], 'bp_method': "minimum_sum", 'ms_scaling_factor': 0.05, 'osd_method': "osd_cs", 'osd_order': 4, 'channel_update': None, 'seed': 42, 'max_iter': 9, 'tqdm_disable' : 1, 'error_bar_precision_cutoff': 1e-6 }

# if len(sys.argv) < 4:

#     print("This script requires inputs l, m, min_k, min_d")

# l = int(sys.argv[1])
# m = int(sys.argv[2])
# min_k = int(sys.argv[3])
# min_d_max = int(sys.argv[4])



# print(f"l = {l}, m = {m}")

# viable_entries = []
# count = 0

# filename = f"bb5_l{l}_m{m}_k{min_k}_d{min_d_max}"

# temp_file = f"../found_codes/{filename}_partial.jsonl"
# results_file = f"../found_codes/{filename}.jsonl"

# # Supprimer ancien fichier temporaire si présent
# if os.path.exists(temp_file):
#     os.remove(temp_file)

# all_ijs = list(itertools.product(range(l), range(m))) # contains every possible (i, j) for this l and m 
# random.seed(42)
# # random.shuffle(all_ijs)

# combos = itertools.product(all_ijs, all_ijs, all_ijs, all_ijs, all_ijs) 

# for combo in combos:

#     ij0, ij1, ij2, ij3, ij4 = combo

#     # Adding the same matrices together would change code stabiliser weight
#     if ij0 == ij1:
#         continue
#     if ij2 == ij3 or ij2 == ij4 or ij3 == ij4:
#         continue

#     Aij = [ij0, ij1]
#     Bij = [ij2, ij3, ij4]
#     # (Could make Aij have three terms and Bij have two but if searching all the terms anyway that will just swap the roles of the left and right data qubits).

#     Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)

#     k = 0

#     if is_valid(Hx, Hz):
#         k = find_k(Hx, Hz)

#     if k < min_k:
#         continue
    
#     p = 0.1
#     target_runs = 500

#     with suppress_stdout():

#         bb5 = css_decode_sim.css_decode_sim(hx = Hx, hz = Hz, error_rate = p, target_runs = target_runs, **osd_options)
#         d_max = bb5.min_logical_weight

#     if d_max < min_d_max:
#         continue

#     entry = {
#         "nkd": [2 * l * m, k, bb5.min_logical_weight],
#         "l" : l,
#         "m" : m,
#         "Aij": Aij,
#         "Bij": Bij,
#         # f"pL at p = {p}": bb5.osdw_logical_error_rate,
#     }

#     print(entry)
    
#     viable_entries.append(entry)
#     count += 1

#     if count % 10 == 0: # save every 10
#         append_entries_to_json(viable_entries, temp_file)
#         viable_entries = []
#         print(f"{count} viable entries saved to {temp_file}")

# # Save any remaining entries after loop
# if viable_entries:
#     append_entries_to_json(viable_entries, temp_file)

# # Rename temp file to final results file
# if os.path.exists(temp_file):
#     os.rename(temp_file, results_file)
#     print(f"Final results saved to {results_file}")



##### GIPPITY code which apparently clears RAM:

# import stim
# import scipy
# import pandas as pd
# import itertools
# import random
# import sys
# import os
# import json
# import gc
# import psutil

# sys.path.append(os.path.abspath("../../src"))
# from bb_ions import *

# def append_entries_to_json(entries, filename):
#     """Append a list of dicts to a JSON lines file."""
#     with open(filename, "a") as f:
#         for entry in entries:
#             entry_clean = {k: list(v) if isinstance(v, tuple) else v for k, v in entry.items()}
#             json.dump(entry_clean, f)
#             f.write("\n")

# # --- Configuration ---
# osd_options = {
#     'xyz_error_bias': [1, 1, 1],
#     'bp_method': "minimum_sum",
#     'ms_scaling_factor': 0.05,
#     'osd_method': "osd_cs",
#     'osd_order': 4,
#     'channel_update': None,
#     'seed': 42,
#     'max_iter': 9,
#     'tqdm_disable': 1,
#     'error_bar_precision_cutoff': 1e-6
# }

# if len(sys.argv) < 5:
#     print("Usage: find_bb5_codes.py l m min_k min_d")
#     sys.exit(1)

# l = int(sys.argv[1])
# m = int(sys.argv[2])
# min_k = int(sys.argv[3])
# min_d_max = int(sys.argv[4])

# print(f"→ Starting (l={l}, m={m}, min_k={min_k}, min_d={min_d_max})")

# # --- Output paths ---
# filename = f"bb5_l{l}_m{m}_k{min_k}_d{min_d_max}"
# temp_file = f"../found_codes/{filename}_partial.jsonl"
# results_file = f"../found_codes/{filename}.jsonl"

# if os.path.exists(temp_file):
#     os.remove(temp_file)

# all_ijs = list(itertools.product(range(l), range(m)))
# # print(all_ijs)

# random.seed(42)
# random.shuffle(all_ijs)

# combos = itertools.product(all_ijs, all_ijs, all_ijs, all_ijs, all_ijs)
# viable_entries = []
# count = 0

# # --- Process monitoring ---
# proc = psutil.Process(os.getpid())

# for combo in combos:
#     # print(combo, end = "\r")
#     ij0, ij1, ij2, ij3, ij4 = combo
#     if ij0 == ij1:
#         continue
#     if ij2 == ij3 or ij2 == ij4 or ij3 == ij4:
#         continue

#     Aij = [ij0, ij1]
#     Bij = [ij2, ij3, ij4]


#     try:
#         Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)
#         if not is_valid(Hx, Hz):
#             continue

#         k = find_k(Hx, Hz)
#         if k < min_k:
#             continue

#         p = 0.1
#         target_runs = 500

#         with suppress_stdout():
#             bb5 = css_decode_sim.css_decode_sim(
#                 hx=Hx, hz=Hz, error_rate=p, target_runs=target_runs, **osd_options
#             )

#         d_max = bb5.min_logical_weight
#         if d_max < min_d_max:
#             continue

#         entry = {
#             "nkd": [2 * l * m, k, d_max],
#             "l": l,
#             "m": m,
#             "Aij": Aij,
#             "Bij": Bij,
#         }

#         print(entry, end = "\r")
        
#         viable_entries.append(entry)
#         count += 1

    
#         append_entries_to_json(viable_entries, temp_file)
#         viable_entries.clear()
#         gc.collect()

#     except Exception as e:
#         print(f"Error on combo {combo}: {e}")

#     # libérer les grands objets C++/NumPy
#     del Hx, Hz
#     if "bb5" in locals():
#         del bb5
#     gc.collect()

# # Sauvegarde finale
# if viable_entries:
#     append_entries_to_json(viable_entries, temp_file)
#     viable_entries.clear()
#     gc.collect()

# if os.path.exists(temp_file):
#     os.rename(temp_file, results_file)
#     print(f"✔ Final results saved to {results_file}")
# else:
#     print("No valid entries found.")



## Gippity code with shuffled i j powers and just sampling 20 million samples (or all the permutations if less that 20 million):

import stim
import scipy
import pandas as pd
import itertools
import random
import sys
import os
import json
import gc
import psutil

sys.path.append(os.path.abspath("../../src"))
from bb_ions import *

def append_entries_to_json(entries, filename):
    """Append a list of dicts to a JSON lines file."""
    with open(filename, "a") as f:
        for entry in entries:
            entry_clean = {k: list(v) if isinstance(v, tuple) else v for k, v in entry.items()}
            json.dump(entry_clean, f)
            f.write("\n")

# --- Configuration ---
osd_options = {
    'xyz_error_bias': [1, 1, 1],
    'bp_method': "minimum_sum",
    'ms_scaling_factor': 0.05,
    'osd_method': "osd_cs",
    'osd_order': 4,
    'channel_update': None,
    'seed': 42,
    'max_iter': 9,
    'tqdm_disable': 1,
    'error_bar_precision_cutoff': 1e-6
}

if len(sys.argv) < 5:
    print("Usage: find_bb5_codes.py l m min_k min_d")
    sys.exit(1)

l = int(sys.argv[1])
m = int(sys.argv[2])
min_k = int(sys.argv[3])
min_d_max = int(sys.argv[4])

print(f"→ Starting (l={l}, m={m}, min_k={min_k}, min_d={min_d_max})")

# --- Output paths ---
filename = f"bb5_l{l}_m{m}_k{min_k}_d{min_d_max}"
temp_file = f"../found_codes/{filename}_partial.jsonl"
results_file = f"../found_codes/{filename}.jsonl"

if os.path.exists(temp_file):
    os.remove(temp_file)

all_ijs = list(itertools.product(range(l), range(m)))
N = len(all_ijs)
T = N**5

# --- Subset sampling parameters ---
desired_samples = 20_000_000  # adjust to fit runtime (~1–2 days)
num_samples = min(desired_samples, T)
seed = 42
random.seed(seed)
max_in_memory = 1_000_000  # threshold to shuffle all combos in memory

# --- Functions ---
def int_to_tuple(x, N):
    return [(x // (N**i)) % N for i in reversed(range(5))]

def sample_random_combos(num_samples, N, T, all_ijs, seed=42):
    """Randomly sample num_samples combos without replacement."""
    seen = set()
    while len(seen) < num_samples:
        idx = random.randrange(T)
        if idx in seen:
            continue
        seen.add(idx)
        yield tuple(all_ijs[i] for i in int_to_tuple(idx, N))

def full_shuffle_combos(all_ijs):
    """Return iterator over all combos shuffled in memory."""
    all_combos = list(itertools.product(all_ijs, repeat=5))
    random.shuffle(all_combos)
    return iter(all_combos)

# --- Choose generation method ---
if T <= max_in_memory:
    print(f"Total combos T={T} small, using full shuffle in memory")
    combos_iter = full_shuffle_combos(all_ijs)
else:
    print(f"Total combos T={T} large, sampling {num_samples} random combos")
    combos_iter = sample_random_combos(num_samples, N, T, all_ijs, seed)

proc = psutil.Process(os.getpid())

for combo in combos_iter:
    ij0, ij1, ij2, ij3, ij4 = combo

    # Early pruning
    if ij0 == ij1:
        continue
    if ij2 == ij3 or ij2 == ij4 or ij3 == ij4:
        continue

    Aij = [ij0, ij1]
    Bij = [ij2, ij3, ij4]

    try:
        Hx, Hz = make_parity_check_matrices(l, m, Aij, Bij)
        if not is_valid(Hx, Hz):
            continue

        k = find_k(Hx, Hz)
        if k < min_k:
            continue

        p = 0.1
        target_runs = 500

        with suppress_stdout():
            bb5 = css_decode_sim.css_decode_sim(
                hx=Hx, hz=Hz, error_rate=p, target_runs=target_runs, **osd_options
            )

        d_max = bb5.min_logical_weight
        if d_max < min_d_max:
            continue

        entry = {
            "nkd": [2 * l * m, k, d_max],
            "l": l,
            "m": m,
            "Aij": Aij,
            "Bij": Bij,
        }

        print(entry, end="\r")
        
        # Save immediately
        append_entries_to_json([entry], temp_file)

    except Exception as e:
        print(f"Error on combo {combo}: {e}")

    # Free large objects
    del Hx, Hz
    if "bb5" in locals():
        del bb5
    gc.collect()

# Final rename
if os.path.exists(temp_file):
    os.rename(temp_file, results_file)
    print(f"✔ Final results saved to {results_file}")
else:
    print("No valid entries found.")



