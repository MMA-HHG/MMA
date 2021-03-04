#!/bin/bash

rm -r * # Clean folder

cp $TESTPATH/*.inp .

rm -r sim_*
rm *.h5

module purge
module load intel/17.2 python/3.6.3
python3 $UNIV_INPUT_PATH/process_multiparametric.py -i-reg ELI5.inp -i-mp ELI5_MP.inp -multiparam-groups -keep-intermediate -univ-inps -ohdf5 results.h5 -g inputs

cp multiparameters/*.h5 .


for simulation in results_*.h5; do

    ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
    mkdir sim_$ksimulation
    mv $simulation sim_$ksimulation/
    cd sim_$ksimulation
        echo "executing simulation $ksimulation"
        $TESTPATH/run_CUPRAD.sh
    cd ..

done;

echo "done" 
