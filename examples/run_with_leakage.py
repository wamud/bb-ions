import deltakit_stim
import sys
sys.modules["stim"] = deltakit_stim

import sinter
import numpy as np
import glob
from stimbposd import SinterDecoder_BPOSD, sinter_decoders
import time
import sys
import os
import multiprocessing


def main():
 
    p = 1e-3
    memory_basis = 'z'
    circuit = deltakit_stim.Circuit.generated(
        f"surface_code:rotated_memory_{memory_basis}",
        rounds=15,
        distance=5,
        after_clifford_depolarization=p)

    tasks = [
        sinter.Task(
            circuit = circuit,
            json_metadata = {
                "p": p,
                "b": memory_basis
            }
        )]
    
    csv_path = f"example_stats/with_leakage.csv"

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 10000,
        max_errors = 100,
        tasks = tasks,
        decoders=['bposd'],
        # existing_data_filepaths = existing,
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
    main()
