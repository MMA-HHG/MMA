#!/bin/bash


test_presets=("test1_modT1" "test1_mod_incT1" "test1_mod_decT1"
              "test1_modT2" "test1_mod_incT2" "test1_mod_decT2"
              "test1_modT3" "test1_mod_incT3" "test1_mod_decT3"
              "test9_modT1" "test9_mod_incT1" "test9_mod_decT1"
              "test9_modT2" "test9_mod_incT2" "test9_mod_decT2"
              "test9_modT3" "test9_mod_incT3" "test9_mod_decT3"
              "test10_modT1" "test10_mod_incT1" "test10_mod_decT1"
              "test10_modT2" "test10_mod_incT2" "test10_mod_decT2"
              "test10_modT3" "test10_mod_incT3" "test10_mod_decT3"
              "test11_modT1" "test11_mod_incT1" "test11_mod_decT1"
              "test11_modT2" "test11_mod_incT2" "test11_mod_decT2"
              "test11_modT3" "test11_mod_incT3" "test11_mod_decT3")



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
