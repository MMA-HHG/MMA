#!/bin/bash
# Submit the first job and save the JobID as JOBONE
JOBPREPROC=$(sbatch --parsable pre_processor.slurm)
# Submit the second job, use JOBONE as depend, save JobID
echo 'preproc submitted'
echo $JOBPREPROC
JOBMAIN=$(sbatch --dependency=afterok:$JOBPREPROC CUPRAD.slurm)
# JOBMAIN=$(sbatch CUPRAD.slurm)
# Submit last job using JOBTWO as depend, do not need to save the JobID

