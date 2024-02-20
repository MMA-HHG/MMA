#!/bin/bash


# test_presets=("test1_modT1" "test1_mod_incT1" "test1_mod_decT1"
#               "test1_modT2" "test1_mod_incT2" "test1_mod_decT2"
#               "test1_modT3" "test1_mod_incT3" "test1_mod_decT3")

test_presets=("test13_modT1" "test13_mod_dec2T1" "test13_mod_dec10T1"
              "test13_modT2" "test13_mod_dec2T2" "test13_mod_dec10T2"
              "test13_modT3" "test13_mod_dec2T3" "test13_mod_dec10T3")



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
