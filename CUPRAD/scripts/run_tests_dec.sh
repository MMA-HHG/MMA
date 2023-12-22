#!/bin/bash



test_presets=("test1_mod_decT1" "test1_mod_decT2" "test1_mod_decT3" "test2_mod_decT1" "test2_mod_decT2" "test2_mod_decT3")



echo ${test_presets[1]}
echo ${test_presets[4]}



module purge
module load hdf5
module load intel-mkl

for k1 in ${!test_presets[@]}; do

echo "----------------------------------------------------------------------"
echo ${test_presets[k1]}
echo "----------------------------------------------------------------------"

mkdir ${test_presets[k1]}
cd ${test_presets[k1]}
$CUPRAD_HOME/build/make_start.e <<INPUTS
${test_presets[k1]}
0
0
0
INPUTS

sbatch $CUPRAD_HOME/scripts/cuprad_JZ32.slurm

cd ..

done


echo "done" 
