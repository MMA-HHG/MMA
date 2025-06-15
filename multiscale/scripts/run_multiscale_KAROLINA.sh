#!/bin/bash

# Submit the pre-processor
JOB1=$(sbatch --parsable $MULTISCALE_SCRIPTS/slurm/KAROLINA/CUPRAD_pre_processor.slurm)

# Submit the main job when the pre-processor is finished
JOB2=$(sbatch --parsable --dependency=afterok:$JOB1 $MULTISCALE_SCRIPTS/slurm/KAROLINA/CUPRAD.slurm)

JOB3=$(sbatch --parsable --dependency=afterok:$JOB2 $MULTISCALE_SCRIPTS/slurm/KAROLINA/CTDSE_prepare_MPI.slurm)

# Run TDSE
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 --nodes=2 $MULTISCALE_SCRIPTS/slurm/KAROLINA/MPI_CTDSE.slurm)

# Collect & merge data
JOB5=$(sbatch --parsable --dependency=afterok:$JOB4 $MULTISCALE_SCRIPTS/slurm/KAROLINA/CTDSE_merge_hdf5.slurm)

# JOB6=$(sbatch --dependency=afterok:$JOB5 $TDSE_1D_SCRIPTS/remove_temp.slurm)


JOBH1=$(sbatch --dependency=afterok:$JOB5 $MULTISCALE_SCRIPTS/slurm/KAROLINA/Hankel_multiprocessing.slurm)