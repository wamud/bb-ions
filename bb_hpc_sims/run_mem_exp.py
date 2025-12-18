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
    
    circuit_paths = glob.glob(f"/288_code_copy/*.stim")

    for path in circuit_paths:
        print(path)

    pbs_jobid = os.environ.get("PBS_JOBID")
    job_number = pbs_jobid.split(".", 1)[0]

    csv_path = f"helios_stats_{job_number}.csv"


    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = 64,
        max_shots = 1_000,
        max_errors = 10,


        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                # # Rebecca's recommended settings:
                # # max_bp_iters = 10, # default 30
                # bp_method="minimum_sum", # product_sum (default), min_sum, min_sum_log
                # ms_scaling_factor = 0.625, # normalisation
                # schedule="serial", 
                # osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                # osd_order=9

                # Long chains settings:
                max_bp_iters = 10_000,
                bp_method = "min_sum",
                osd_order = 5,
                osd_method = "osd_cs"
                # (other settings left unspecified so they take the default values)
                
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
