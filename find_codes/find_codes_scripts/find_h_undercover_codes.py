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


def factor_pairs(n):
    """Returns all factor pairs (u,t) with u*t=n."""
    pairs = []
    for i in range(1, int(n**0.5) + 1):
        if n % i == 0:
            pairs.append((i, n // i))
    return pairs


def append_entries_to_json(entries, filename):
    with open(filename, "a") as f:
        for entry in entries:
            entry_clean = {
                k: list(v) if isinstance(v, tuple) else v
                for k, v in entry.items()
            }
            json.dump(entry_clean, f)
            f.write("\n")


def undercover(Aij, Bij, l, m, u, t):
    """
    Given a code (Aij,Bij,l,m), assume it is a (u*t)-cover.
    Return the unique candidate base code.
    """

    base_l = l // u
    base_m = m // t

    A_base = [(i % base_l, j % base_m) for i, j in Aij]
    B_base = [(i % base_l, j % base_m) for i, j in Bij]

    return A_base, B_base, base_l, base_m


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
    print(
        "Usage:\n"
        "python find_h_undercover_codes.py gross_code h min_k min_d"
    )
    sys.exit()

code_name = str(sys.argv[1])
h = int(sys.argv[2])
min_k = int(sys.argv[3])
min_d = int(sys.argv[4])

code = eval(f"{code_name}()")

print(
    f"Finding possible base codes for which "
    f"[[{code.n}, {code.k}, {code.d_max}]] "
    f"is a {h}-cover."
)

og_l = code.l
og_m = code.m
og_Aij = code.Aij
og_Bij = code.Bij

filename = (
    f"[[{code.n}, {code.k}, {code.d_max}]] "
    f"{h}-undercovers "
    f"min_k={min_k}, min_d={min_d}"
)

temp_file = f"../found_h_undercover_codes/{filename}_partial.jsonl"
results_file = f"../found_h_undercover_codes/{filename}.jsonl"

pairs = factor_pairs(h)
pairs = list(set(pairs + [(b, a) for a, b in pairs]))

seen_structures = set()

for u, t in pairs:
    print(u,t)

    # Base lattice must divide the current lattice.
    if og_l % u != 0:
        continue
    if og_m % t != 0:
        continue

    Aij, Bij, l, m = undercover(
        og_Aij,
        og_Bij,
        og_l,
        og_m,
        u,
        t,
    )

    # Remove translation symmetry
    imin, jmin = Aij[0]

    Aij_norm = [
        ((i - imin) % l, (j - jmin) % m)
        for i, j in Aij
    ]

    Bij_norm = [
        ((i - imin) % l, (j - jmin) % m)
        for i, j in Bij
    ]

    Aij_norm = sorted(Aij_norm)
    Bij_norm = sorted(Bij_norm)

    key = (tuple(Aij_norm), tuple(Bij_norm))

    if key in seen_structures:
        continue

    seen_structures.add(key)

    # Duplicate terms mean this cannot be a valid base code.
    if len(set(Aij_norm)) != len(Aij_norm):
        continue

    if len(set(Bij_norm)) != len(Bij_norm):
        continue

    Hx, Hz = make_parity_check_matrices(
        l,
        m,
        Aij_norm,
        Bij_norm
    )

    if not is_valid(Hx, Hz):
        continue

    k = find_k(Hx, Hz)

    if k < min_k:
        continue

    p = 0.1
    target_runs = 1000

    with suppress_stdout():
        bb = css_decode_sim.css_decode_sim(
            hx=Hx,
            hz=Hz,
            error_rate=p,
            target_runs=target_runs,
            **osd_options
        )

        d_max = bb.min_logical_weight

    if d_max < min_d:
        del bb
        gc.collect()
        continue

    entry = {
        "nkd": [2 * l * m, k, d_max],
        "l": l,
        "m": m,
        "Aij": Aij_norm,
        "Bij": Bij_norm,
        "cover_factorisation": {
            "u": u,
            "t": t,
            "h": h
        }
    }

    print(entry)

    append_entries_to_json([entry], temp_file)

    del bb
    gc.collect()


if os.path.exists(temp_file):
    os.rename(temp_file, results_file)
    print(f"Final results saved to {results_file}")