#!/bin/bash



test_presets=("test3_modT1" "test3_mod_incT1" "test3_mod_decT1"
              "test3_modT2" "test3_mod_incT2" "test3_mod_decT2"
              "test3_modT3" "test3_mod_incT3" "test3_mod_decT3"
              "test4_modT1" "test4_mod_incT1" "test4_mod_decT1"
              "test4_modT2" "test4_mod_incT2" "test4_mod_decT2"
              "test4_modT3" "test4_mod_incT3" "test4_mod_decT3")



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
