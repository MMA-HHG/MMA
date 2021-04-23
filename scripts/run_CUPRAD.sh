#!/bin/bash

# Submit the pre-processor
JOBPREPROC=$(sbatch --parsable $TESTPATH/pre_processor.slurm)

# Submit the main job when the pre-processor is finished
JOBMAIN=$(sbatch --dependency=afterok:$JOBPREPROC $TESTPATH/CUPRAD.slurm)
