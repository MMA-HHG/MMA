#!/bin/bash

# Prepare starting files

# JOBPREPARE=$(qsub -N preparejobs submit_prepare)

rm -r sim_*
rm *.h5

var=*.slurm
echo $var

module purge
module load intel/17.2 python/3.6.3

python3 $UNIV_INPUT/process_multiparametric.py -i-reg ELI1.inp -i-mp ELI1_MP.inp -univ-inps -ohdf5 reults.h5 -g inputs

cp multiparametric/*.h5 .


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
