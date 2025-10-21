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
    
    empty_folder("../collected_stats") ## CURRENTLY DELETING COLLECTED_STATS EACH TIME TO TEST TIME IT TAKES GIVEN DIFFERENT NUMBER OF CORES USED. REMOVE THIS LINE WHEN ACTUALLY RUNNING! DON"T WANNA DELETE ALL YOUR STATS
    start_time = time.time()
    
    circuit_paths = glob.glob(f"../circuits/*1000*.stim")
    csv_path = f"../collected_stats/collected_stats.csv"

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = 4,
        max_shots = 4,
        max_errors = 4,
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
