#!/bin/bash

# Submit the pre-processor
JOBPREPROC=$(sbatch --parsable $CUPRAD_SCRIPTS/pre_processor.slurm)

# Submit the main job when the pre-processor is finished
JOBMAIN=$(sbatch --dependency=afterok:$JOBPREPROC $CUPRAD_SCRIPTS/CUPRAD.slurm)

# Create folder for TDSE

# Run TDSE

# Collect & merge data

# Run Hankels


