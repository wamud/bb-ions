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
        print("Usage: run_mem_exp.py <num_workers> <suffix>")
        sys.exit(1)

    num_workers = int(sys.argv[1])
    csv_file_number = sys.argv[2]    
 
    start_time = time.time()
    
    circuit_paths = glob.glob(f"../circuits/*12*.stim")
    csv_path = f"../collected_stats/collected_stats_{csv_file_number}.csv"

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = num_workers,
        max_shots = 5000,
        max_errors = 5000,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                # max_bp_iters = 10,
                bp_method="minimum_sum", # product_sum, min_sum, min_sum_log
                ms_scaling_factor = 0.62, # normalisation
                schedule="serial",
                osd_method="osd0", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD
                osd_order=0
            )
        },
        print_progress = True
        )

    end_time = time.time()
    print(f"Finished collecting in {(end_time - start_time):.2f} seconds")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
