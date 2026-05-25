import deltakit_stim
# injections = ['stim._detect_machine_architecture', 'stim._stim_polyfill', 'stim']
import sys
# for namespace in injections:
#     sys.modules[namespace] = sys.modules[f"deltakit_{namespace}"] # setting the stim one equal to the deltakit_stim one


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

 
    start_time = time.time()

    print(f"Start time = {start_time}")
    

    circuit_paths = glob.glob(f"../circuits/with_leakage/helios/*.stim") 
    circuit_paths.sort()

    for path in circuit_paths:
        print(path)

    if len(circuit_paths) == 0:
        print("No circuits")
    
    csv_path = f"../collected_stats/helios_noise/other_investigations/gross_w_leakage.csv"

    # existing = glob.glob(f"../collected_stats/helios_noise/other_investigations/*parallel*.stim") 

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]


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

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 100_000,
        max_errors = 10,
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

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 1_000_000,
        max_errors = 2,
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

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 10_000_000,
        max_errors = 2,
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


    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 100_000_000,
        max_errors = 2,
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

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 1_000_000_000,
        max_errors =2,
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

    end_time = time.time()
    print(f"Finished collecting in {(end_time - start_time):.2f} seconds")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
