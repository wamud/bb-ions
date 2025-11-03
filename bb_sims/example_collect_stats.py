#!/usr/bin/python3

import stim
import pymatching
import sinter
from typing import *


def main():
    surface_code_tasks = [
        sinter.Task(
            circuit = stim.Circuit.generated(
                f"surface_code:rotated_memory_{b}",
                rounds=d * 3,
                distance=d,
                after_clifford_depolarization=noise,
                after_reset_flip_probability=noise,
                before_measure_flip_probability=noise,
                before_round_data_depolarization=noise,
            ),
            json_metadata={'d': d, 'r': d * 3, 'p': noise,'mem':b},
            )
        for d in [3]
        for noise in [0.001]
        for b in 'z'
        ]

    collected_surface_code_stats: List[sinter.TaskStats] = sinter.collect(
            num_workers = 1,   # (num cores)
        tasks=surface_code_tasks,
        decoders=['pymatching'],
        max_shots=100,
        max_errors=100,
        print_progress=False,
        save_resume_filepath='collected_stats.csv',
    )

if __name__ == "__main__":
    import multiprocessing
    main()
