import deltakit_stim
import sys
# sys.modules["stim"] = deltakit_stim 

import sinter
import numpy as np
import glob
from stimbposd import SinterDecoder_BPOSD, sinter_decoders
import time
import sys
import os
import multiprocessing


def main():

    circuit = deltakit_stim.Circuit("""
    R 0 1
    LEAKAGE(0.1) 0 1
    CNOT 0 1
    HERALDED_ERASE(0.001) 0 1
    M 0 1
    DETECTOR rec[-1]
    DETECTOR rec[-2]
    HERALD_LEAKAGE_EVENT 0 1
    DETECTOR rec[-1]
    DETECTOR rec[-2]
    OBSERVABLE_INCLUDE(0) rec[-1]
    """)
  
    csv_path = f"example_stats/leaky_lossy.csv"

    tasks = [sinter.Task(circuit=circuit, json_metadata={"p": 0.1, "circ": "example"})]

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 100000,
        max_errors = 10,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                max_bp_iters = 10_000, # default 30
                bp_method = "min_sum", # product_sum (default), min_sum, min_sum_log
                # ms_scaling_factor = 0.625, # normalisation
                # schedule = "serial", 
                osd_method = "osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                osd_order = 5 
            )
        },
        print_progress = True
        )



if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
