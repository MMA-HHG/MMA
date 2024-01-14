#!/bin/bash


test_presets=("test14_mod+100" "test14_mod+80" "test14_mod+60"
              "test14_mod+40" "test14_mod+20" "test14_mod0"
              "test14_mod-20" "test14_mod-40" "test14_mod-60"
              "test14_mod-80" "test14_mod-100")



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
