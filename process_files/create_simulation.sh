#!/bin/bash

#module purge
#module load compiler/intel

cd sims

n=0
while ! mkdir test$n
do
	n=$((n+1))
done
cd ..

cp FreeFormInputs.inp sims/test$n
cp run_occigen2.slurm sims/test$n


cd sims/test$n

sbatch run_occigen2.slurm -J test$n

