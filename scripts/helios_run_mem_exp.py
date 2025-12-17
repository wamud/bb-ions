import stim
import sinter
import numpy as np
import glob
from stimbposd import SinterDecoder_BPOSD, sinter_decoders
import time
import sys
import os

def main():

 
    start_time = time.time()
    

    circuit_paths = glob.glob(f"../circuits/helios/exclude_opp_basis_detectors/*.stim")
    circuit_paths.sort()

    
    csv_path = f"../collected_stats/helios_noise_long_chain_BPOSD_settings.csv"
    # existing = [f"../collected_stats/tham_modules_noise_long_chain_BPOSD_settings_incl_opp_detectors.csv"]

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]


    samples = sinter.collect(
        num_workers = 64,
        max_shots = 1000,
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
        num_workers = 64,
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
        print_progress = False
        )

    samples = sinter.collect(
        num_workers = 64,
        max_shots = 100000,
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
        print_progress = False
        )

    samples = sinter.collect(
        num_workers = 64,
        max_shots = 1_000_000,
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
        print_progress = False
        )

    samples = sinter.collect(
        num_workers = 64,
        max_shots = 10_000_000,
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
        print_progress = False
        )


    samples = sinter.collect(
        num_workers = 64,
        max_shots = 1_000_000_000_000,
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
        print_progress = False
        )

    end_time = time.time()
    print(f"Finished collecting in {(end_time - start_time):.2f} seconds")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
