#!/bin/bash



test_presets=("test1_modT1" "test1_mod_incT1" "test1_mod_decT1"
              "test1_modT2" "test1_mod_incT2" "test1_mod_decT2"
              "test1_modT3" "test1_mod_incT3" "test1_mod_decT3"
              "test5_modT1" "test5_mod_incT1" "test5_mod_decT1"
              "test5_modT2" "test5_mod_incT2" "test5_mod_decT2"
              "test5_modT3" "test5_mod_incT3" "test5_mod_decT3"
              "test6_modT1" "test6_mod_incT1" "test6_mod_decT1"
              "test6_modT2" "test6_mod_incT2" "test6_mod_decT2"
              "test6_modT3" "test6_mod_incT3" "test6_mod_decT3"
              "test7_modT1" "test7_mod_incT1" "test7_mod_decT1"
              "test7_modT2" "test7_mod_incT2" "test7_mod_decT2"
              "test7_modT3" "test7_mod_incT3" "test7_mod_decT3"
              "test8_modT1" "test8_mod_incT1" "test8_mod_decT1"
              "test8_modT2" "test8_mod_incT2" "test8_mod_decT2"
              "test8_modT3" "test8_mod_incT3" "test8_mod_decT3"
              "test9_modT1" "test9_mod_incT1" "test9_mod_decT1"
              "test9_modT2" "test9_mod_incT2" "test9_mod_decT2"
              "test9_modT3" "test9_mod_incT3" "test9_mod_decT3")


# test_presets=("test8_modT1" "test8_mod_incT1" "test8_mod_decT1"
#               "test8_modT2" "test8_mod_incT2" "test8_mod_decT2"
#               "test8_modT3" "test8_mod_incT3" "test8_mod_decT3")



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
