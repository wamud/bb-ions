#!/bin/bash

# Run this as: qsub submit.sh

#PBS -N collect_stats
#PBS -l ncpus=1
#PBS -l mem=10GB
#PBS -l walltime=00:01:00
#PBS -m abe
#PBS -M anthony.orourke@student.uts.edu.au

###################################
# Setup any input files for the run
###################################


# Create a /scratch directory with a unique name.
SCRATCH="/scratch/${USER}_${PBS_JOBID%.*}"
mkdir -p ${SCRATCH}

cd ${PBS_O_WORKDIR}

# Copy your input files from there to the scratch directory you created above.
if [ -f collected_stats.csv ]; then
    cp collected_stats.csv ${SCRATCH}
fi

###############
# Start the Job
###############

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bb_env

cd ${SCRATCH}
python3 ${PBS_O_WORKDIR}/run_mem_exp.py $1 

#####################################################
# Copy results back to your own directory and cleanup
#####################################################

mv ${SCRATCH}/collected_stats.csv ${PBS_O_WORKDIR}

cd ${PBS_O_WORKDIR}
rmdir ${SCRATCH}

