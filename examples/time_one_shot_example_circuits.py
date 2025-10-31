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
    csv_suffix  = sys.argv[2]    

    for num_rounds in [60, 80, 90, 100]:

        print(num_rounds)
        
        empty_folder("example_stats")

        csv_path = f"example_stats/collected_stats_{csv_suffix}.csv"

        circuit_paths = glob.glob(f"example_circuits/*r={num_rounds},*.stim")        
        

        start_time = time.time()


        tasks = [
            sinter.Task(
                circuit_path = path,
                json_metadata = sinter.comma_separated_key_values(path),
            )
            for path in circuit_paths
        ]

        samples = sinter.collect(
            num_workers = num_workers,
            max_shots = 1,  
            max_errors = 1,
            tasks = tasks,
            decoders=['bposd'],
            save_resume_filepath = csv_path,
            custom_decoders = {
                "bposd": SinterDecoder_BPOSD(
                    # max_bp_iters = 10, # default is 30
                    bp_method="minimum_sum", # product_sum (default), min_sum, min_sum_log
                    ms_scaling_factor = 0.625, # normalisation
                    schedule="serial",
                    osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                    osd_order=9 # default is 60
                )
            },
            print_progress = True
            )

        end_time = time.time()
        print(f"  Finished in {(end_time - start_time):.2f} seconds\n")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
