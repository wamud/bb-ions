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

sys.path.append(os.path.abspath("../packages/BeamSearchDecoder"))
from beamsearch import BeamSearch
from sinter_beamsearch import SinterDecoder_BeamSearch



custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                ## Original BB code paper (and subsequent long chains paper) settings:
                max_bp_iters = 10_000, # default 30
                bp_method = "min_sum", # product_sum (default), min_sum, min_sum_log
                osd_method = "osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                osd_order = 5,
                # Implies ms_scaling_factor is 1.0, schedule is parallel
            ),
            "bposd_625": SinterDecoder_BPOSD(
                ## Rebecca's settings (slower but slightly better)
                max_bp_iters = 10, # Becca: 10, # default 30
                bp_method="min_sum", # product_sum (default), min_sum, min_sum_log
                ms_scaling_factor = 0.625, # normalisation
                schedule="serial", 
                osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                osd_order=9
            ),
            
            # Beamsearch decoders:
            "beam8_230iters": SinterDecoder_BeamSearch(),
            "beam32_340iters": SinterDecoder_BeamSearch(beam_width=32, initial_iters=40, iters_per_round=30),
            "beam64_640iters": SinterDecoder_BeamSearch(max_rounds=20, beam_width=64, initial_iters=40, iters_per_round=30),
            "beam64_32res_640iters": SinterDecoder_BeamSearch(
                max_rounds=20, beam_width=64, num_results=32, initial_iters=40, iters_per_round=30
            ),
        }

def main():

 
    

    circuit_paths = glob.glob(f"../circuits/leakage_and_loss/all_codes/*.stim") 
    csv_path = f"../collected_stats/helios_noise/leakage_and_loss/testing_stats.csv"

    circuit_paths.sort()
    if len(circuit_paths) == 0:
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!No circuits!!!!!!!!!!!!!!!!!!!!!!!!!!")
    for path in circuit_paths:
        print(os.path.basename(path))

    existing =  [] #glob.glob(f"../collected_stats/helios_noise/leakage_and_loss/*new_all_codes*.csv") 
    # print(f"existing = {existing}")

    tasks = [
        sinter.Task(
            circuit = deltakit_stim.Circuit.from_file(path),#.flattened(), # Have to flatten to make leakage heralds work, otherwise for some reason (deltakit stim bug?) 4 or more rounds in a circuit create a mismatch between the num. of detectors in the circuit versus the DEM.
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    start_time = time.time()
    print(f"Start time = {start_time}")

    samples = sinter.collect(
        num_workers = multiprocessing.cpu_count(),
        max_shots = 1000,
        max_errors = 5,
        tasks = tasks,
        decoders=['beam64_32res_640iters'],
        existing_data_filepaths = existing,
        save_resume_filepath = csv_path,
        custom_decoders = custom_decoders,
        print_progress = True
        )


    end_time = time.time()
    print(f"Finished collecting in {(end_time - start_time):.2f} seconds")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
