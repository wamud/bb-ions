#!/bin/bash

# Run this as: qsub submit.sh

#PBS -N sims
#PBS -l ncpus=64
#PBS -l mem=300GB
#PBS -l walltime=336:00:00
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
if [ -d "../circuits/helios/exclude_opp_basis_detectors/288_code_copy/" ]; then
    cp -r "../circuits/helios/exclude_opp_basis_detectors/288_code_copy/" "${SCRATCH}"
fi


###############
# Start the Job
###############

source ~/miniconda3/etc/profile.d/conda.sh
conda activate bb_env

cd ${SCRATCH}
python3 ${PBS_O_WORKDIR}/run_mem_exp.py

#####################################################
# Copy results back to your own directory and cleanup
#####################################################

mv ${SCRATCH}/collected_stats*.csv ${PBS_O_WORKDIR}

cd ${PBS_O_WORKDIR}
rm -rf ${SCRATCH}

