#!/bin/bash

# Submit the pre-processor
JOBONE=$(sbatch --parsable $CUPRAD_SCRIPTS/pre_processor.slurm)
touch 1.test

# Submit the main job when the pre-processor is finished
JOBTWO=$(sbatch --dependency=afterok:$JOBONE $CUPRAD_SCRIPTS/CUPRAD.slurm)

# Create folder for TDSE & prepare input
mkdir TDSEs
cd TDSEs
touch 2.test

# Run TDSE

# Collect & merge data

# Run Hankels


