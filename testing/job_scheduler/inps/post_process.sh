#!/bin/bash

rm -r * # Clean folder

sbatch $TESTPATH/plot1.slurm

echo "done" 
