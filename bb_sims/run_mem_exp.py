import stim
import sinter
import numpy as np
import glob
from stimbposd import SinterDecoder_BPOSD, sinter_decoders
import time
import sys



def main():
    if len(sys.argv) < 3:
        print("Usage: run_mem_exp.py <n_k_d> <T2>")
        sys.exit(1)

    nkd = sys.argv[1]
    T2 = int(sys.argv[2])
 
    start_time = time.time()
    
    # circuit_paths = glob.glob(f"{nkd} w T2 = {T2}/pause_0/*.stim")
    

    # # Excluding p=0.0005 circuits:
    circuit_paths = [
        path for path in glob.glob(f"{nkd} w T2 = 10/*/*.stim")
        if "p=0.0005" not in path
    ]


    if len(circuit_paths) == 0:
        print("No circuits")
        sys.exit(1)
        
    csv_path = f"collected_stats_{nkd}__T2={T2}.csv"

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = 8,
        #max_shots = 40_000_000,
        #max_errors = 100,
        max_shots = 1,
        max_errors = 1,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = {
            "bposd": SinterDecoder_BPOSD(
                # max_bp_iters = 10, # default 30
                bp_method="minimum_sum", # product_sum (default), min_sum, min_sum_log
                ms_scaling_factor = 0.625, # normalisation
                schedule="serial", 
                osd_method="osd_cs", # "osd0" - zero-order OSD, "osd_e" - exhaustive OSD, "osd_cs": combination-sweep OSD (default)
                osd_order=9
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
