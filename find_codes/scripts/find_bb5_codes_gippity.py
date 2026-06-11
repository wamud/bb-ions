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


filename = f"bb5_l{l}_m{m}_k{min_k}_d{min_d_max}"

temp_file = f"../found_codes/{filename}_partial.jsonl"
results_file = f"../found_codes/{filename}.jsonl"
progress_file = f"../found_codes/l{l}_m{m}_progress.txt"



# We will make ivectors containing i0, i1, i2, i3, i4 and jvectors containing j0, j1, j2, j3, j4
# which are the powers of the terms in the matrices A and B, i.e. 
# A = x^i0 * y^j0 + x^i1 * y^j1
# B = x^i2 * y^j2 + x^i3 * y^j3 + x^i4 * y^j4
# In turn, A and B will be used in the Bicycle Bivariate code's parity check matrices as Hx = [A|B] and Hz = [B^T|A^T]

random.seed(42)

ivalues = range(l)
ivectors = list(itertools.product(ivalues, repeat = 5))
# random.shuffle(ivectors)

jvalues = range(m)
jvectors = list(itertools.product(jvalues, repeat = 5))
# random.shuffle(jvectors)

# Begin gippity:

from math import gcd

def units_mod(n):
    """Return all invertible elements of Z_n."""
    return [a for a in range(n) if gcd(a, n) == 1]


def canonical_bb(Aij, Bij, l, m):
    """
    Return a canonical representative of a BB code under:

      - common translation
      - automorphisms of Z_l × Z_m

    This avoids testing codes that differ only by coordinate relabelling.
    """

    reps = []

    units_l = units_mod(l)
    units_m = units_mod(m)

    for a in units_l:
        for b in units_m:

            # Apply automorphism:
            # (i,j) -> (a*i mod l, b*j mod m)
            A = sorted(
                ((a * i) % l, (b * j) % m)
                for i, j in Aij
            )

            B = sorted(
                ((a * i) % l, (b * j) % m)
                for i, j in Bij
            )

            # Remove global translation
            ai0, aj0 = A[0]

            A_norm = sorted(
                (
                    (i - ai0) % l,
                    (j - aj0) % m
                )
                for i, j in A
            )

            B_norm = sorted(
                (
                    (i - ai0) % l,
                    (j - aj0) % m
                )
                for i, j in B
            )

            reps.append(
                (
                    tuple(A_norm),
                    tuple(B_norm)
                )
            )

    return min(reps)


# ------------------------------------------------------------------
# Existing search loop
# ------------------------------------------------------------------

seen_structures = set()

for loop in range(len(ivectors)):
    for count, jvec in enumerate(jvectors):
        ivec = ivectors[(count + loop) % len(ivectors)]

        if count % 5000 == 0:
            with open(progress_file, "w") as f:
                f.write(f"loop = {loop}, count = {count}")

        i0 = ivec[0]
        i1 = ivec[1]
        i2 = ivec[2]
        i3 = ivec[3]
        i4 = ivec[4]

        j0 = jvec[0]
        j1 = jvec[1]
        j2 = jvec[2]
        j3 = jvec[3]
        j4 = jvec[4]

        # (Could make Aij have three terms and Bij have two but if searching all the terms anyway that will just swap the roles of the left and right data qubits).
        # Skip values where the same term appears twice in the same matrix (this would change the weight of the stabilisers as when they're added together mod 2 their ones cancel out in the parity check matrices)
        if (i0, j0) == (i1, j1):
            continue

        if (i2, j2) == (i3, j3) or (i2, j2) == (i4, j4) or (i3, j3) == (i4, j4):
            continue

        # Also equivalent polynomials with reordered terms need not be repeated. I.e. x^1y^2 + x^3 is equivalent to x^3 + x^1y^2.
        # To avoid repeats let's just search terms that are in ascending order lexicographically (compare first element of tuple, if equal then compare second element of tuple)
        if (i0, j0) > (i1, j1):
            continue

        if not ((i2, j2) < (i3, j3) < (i4, j4)):
            continue

        Aij = [(i0, j0), (i1, j1)]
        Bij = [(i2, j2), (i3, j3), (i4, j4)]

        # Also quotient by automorphisms of Z_l × Z_m
        # and by common translations of all terms.
        #
        # This avoids testing codes whose parity check
        # matrices differ only by a relabelling of coordinates.

        key_struct = canonical_bb(Aij, Bij, l, m)

        # --- skip structures déjà vues ---
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
            bb5 = css_decode_sim.css_decode_sim(
                hx=Hx,
                hz=Hz,
                error_rate=p,
                target_runs=target_runs,
                **osd_options
            )

            d_max = bb5.min_logical_weight

        if d_max < min_d_max:
            del bb5
            gc.collect()
            continue

        entry = {
            "nkd": [2 * l * m, k, bb5.min_logical_weight],
            "l": l,
            "m": m,
            "Aij": Aij,
            "Bij": Bij,
        }

        print(entry)

        append_entries_to_json([entry], temp_file)

        del bb5
        gc.collect()

# End gippity

# Rename temp file to final results file
if os.path.exists(temp_file):
    os.rename(temp_file, results_file)
    print(f"Final results saved to {results_file}")

