import stim
import sinter
import numpy as np
import glob
from stimbposd import SinterDecoder_BPOSD, sinter_decoders
import time
import sys
import os
sys.path.append(os.path.abspath("../src"))
from bb_ions import *


def main():
    if len(sys.argv) < 2:
        print("Usage: run_mem_exp.py <codeT2>")
        sys.exit(1)

    code_T2 = sys.argv[1]
 
    start_time = time.time()
    
    circuit_paths = glob.glob(f"../circuits/uniform_plus_shift_and_shuttle_w_dephasing_idling/{code_T2}/*.stim")
    csv_path = f"../collected_stats/collected_stats_{code_T2}.csv"

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = 64,
        #max_shots = 40_000_000,
        #max_errors = 100,
        max_shots = 100,
        max_errors = 100,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                # max_bp_iters = 10, # default 30
                bp_method="minimum_sum", # product_sum (default), min_sum, min_sum_log
                ms_scaling_factor = 0.625, # normalisation
                schedule="serial", 
                osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                osd_order=9
            )
        },
        print_progress = False
        )

    end_time = time.time()
    print(f"Finished collecting in {(end_time - start_time):.2f} seconds")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
