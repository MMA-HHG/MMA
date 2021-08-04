#!/bin/bash


for fold in sim_*; do

    # ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
    cd $fold
        echo $fold
        cp results_*.h5 ../
        cp CUPRAD-*.output ../
        cp pre_processor-*.output ../
    cd ..

done;



echo "done" 
