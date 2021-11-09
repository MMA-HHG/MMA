#!/bin/bash

# Submit the pre-processor
JOB1=$(sbatch --parsable $CUPRAD_SCRIPTS/pre_processor.slurm)
# touch 1.test

# Submit the main job when the pre-processor is finished
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 $CUPRAD_SCRIPTS/CUPRAD.slurm)

# Create folder for TDSE & prepare input
mkdir TDSEs
cd TDSEs
# touch 2.test
JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 $TDSE_1D_SCRIPTS/prepare_TDSE.slurm)

# Run TDSE
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 $TDSE_1D_HOME/slurm/scale3_v5-24h.slurm)

# Collect & merge data
mkdir temp
JOB5=$(sbatch --parsable --dependency=afterok:$JOB4 $TDSE_1D_SCRIPTS/merge_all_move.slurm)

JOB6=$(sbatch --dependency=afterok:$JOB5 $TDSE_1D_SCRIPTS/remove_temp.slurm)

# Run Hankels
python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $HANKEL_HOME/process_files/Hankel_Inputs_mp_large_all_cummulative.inp -ohdf5 inputs_Hankel_all_cummulative.h5 -g inputs
python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $HANKEL_HOME/process_files/Hankel_Inputs_mp_large_all_cummulative_noabs.inp -ohdf5 inputs_Hankel_all_cummulative_noabs.h5 -g inputs
python3 $UNIV_INPUT_PATH/create_universal_HDF5.py -i $HANKEL_HOME/process_files/Hankel_Inputs_mp_large_all_cummulative_nodispersion.inp -ohdf5 inputs_Hankel_all_cummulative_nodispersion.h5 -g inputs

JOBH1=$(sbatch --dependency=afterok:$JOB5 $HANKEL_HOME/process_files/Hankel_mp_all_cummulative.slurm)
JOBH2=$(sbatch --dependency=afterok:$JOB5 $HANKEL_HOME/process_files/Hankel_mp_all_cummulative_noabs.slurm)
JOBH3=$(sbatch --dependency=afterok:$JOB5 $HANKEL_HOME/process_files/Hankel_mp_all_cummulative_nodispersion.slurm)

