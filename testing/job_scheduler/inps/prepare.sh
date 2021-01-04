#!/bin/bash

# Prepare starting files

# JOBPREPARE=$(qsub -N preparejobs submit_prepare)

rm -r sim_*

var=*.slurm
echo $var

for fold in results_*.h5; do 

echo $fold
echo $fold | grep -Po '(?<=_)\d+' 
var=$(echo $fold | grep -Po '(?<=_)\d+')
echo $var

mkdir sim_$var
cp $fold sim_$var/
# cp pre_processor.sh sim_$var/
cp pre_processor.slurm sim_$var/
cp CUPRAD.slurm sim_$var/
cp run_CUPRAD.sh sim_$var/

cd sim_$var
echo "runscript"
# ./../pre_processor.sh
./pre_processor.slurm
cd ..


# cd $fold/

# export fold JOBPREPARE;
# ./run_code
# # qsub -N ArPropagationIR_$fold submit_avakas

# echo $fold

# cd ..

done;



echo "done" 
