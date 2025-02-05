#!/bin/bash

# Submit the first job
JOB1=$(qsub $MULTISCALE_SCRIPTS/pbs/CUPRAD_pre_processor.pbs)

# Submit the main job when the first job is finished (depend=afterok)
JOB2=$(qsub -W depend=afterok:$JOB1 $MULTISCALE_SCRIPTS/pbs/CUPRAD.pbs)

JOB3=$(qsub -W depend=afterok:$JOB2 $MULTISCALE_SCRIPTS/pbs/CTDSE_prepare_MPI.pbs)

# Include whatever resources you need either in the PBS script itself or via command-line options
JOB4=$(qsub -W depend=afterok:$JOB3 \
            -l select=1:ncpus=1800 \
            $MULTISCALE_SCRIPTS/pbs/MPI_CTDSE.pbs)

JOB5=$(qsub -W depend=afterok:$JOB4 $MULTISCALE_SCRIPTS/pbs/CTDSE_merge_hdf5.pbs)

# Run Hankel code after the merge completes
JOBH1=$(qsub -W depend=afterok:$JOB5 $MULTISCALE_SCRIPTS/pbs/Hankel_multiprocessing.pbs)