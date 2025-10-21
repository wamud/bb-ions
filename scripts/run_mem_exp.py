import stim
import sinter
import numpy as np
import glob
from stimbposd import SinterDecoder_BPOSD, sinter_decoders


def main():
    circuit_paths = glob.glob(f"../circuits/*.stim")
    csv_path = f"../collected_stats/collected_100_stats.csv"

    tasks = [
        sinter.Task(
            circuit_path = path,
            json_metadata = sinter.comma_separated_key_values(path),
        )
        for path in circuit_paths
    ]

    samples = sinter.collect(
        num_workers = 1,
        max_shots = 1,
        max_errors = 1,
        tasks = tasks,
        decoders=['bposd'],
        save_resume_filepath = csv_path,
        custom_decoders = sinter_decoders(),
        print_progress = True
        )

    print("Collection terminé")


if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()   # utile sous Windows/macOS
    main()
