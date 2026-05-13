import multiprocessing
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

 
    start_time = time.time()
    
    # circuit_paths = glob.glob(f"../circuits/uniform_plus_shift_and_shuttle_w_dephasing_idling/*T2 = 10*/pause_0/*.stim")
    
    # # Excluding 288 code and p=0.0005 circuits:
    # circuit_paths = [
    #     path for path in glob.glob("../circuits/uniform_plus_shift_and_shuttle_w_dephasing_idling/*T2 = 10*/pause_0/*.stim")
    #     if "288_12_18" not in path and "p=0.0005" not in path
    # ]

    circuit_paths = glob.glob(f"../circuits/uniform/*.stim")

    csv_path = f"../collected_stats/uniform/uniform.csv"

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 1000,
        max_errors = 100,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                ## Rebecca's settings:
                # max_bp_iters = 10, # default 30
                # bp_method="minimum_sum", # product_sum (default), min_sum, min_sum_log
                # ms_scaling_factor = 0.625, # normalisation
                # schedule="serial", 
                # osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                # osd_order=9
                max_bp_iters = 10_000, 
                bp_method = "min_sum", 
                osd_method = "osd_cs", 
                osd_order = 5 
            )
        },
        print_progress = True
        )

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 10000,
        max_errors = 100,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                max_bp_iters = 10_000, 
                bp_method = "min_sum", 
                osd_method = "osd_cs", 
                osd_order = 5 
            )
        },
        print_progress = True
        )

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 10000,
        max_errors = 100,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                max_bp_iters = 10_000, 
                bp_method = "min_sum", 
                osd_method = "osd_cs", 
                osd_order = 5 
            )
        },
        print_progress = True
        )
    print("Now running to see a million shots or two errors...")
    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 10000000,
        max_errors = 5,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                # max_bp_iters = 10, # default 30
                # bp_method="minimum_sum", # product_sum (default), min_sum, min_sum_log
                # ms_scaling_factor = 0.625, # normalisation
                # schedule="serial", 
                # osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                # osd_order=9
                max_bp_iters = 10_000, 
                bp_method = "min_sum", 
                osd_method = "osd_cs", 
                osd_order = 5 
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
