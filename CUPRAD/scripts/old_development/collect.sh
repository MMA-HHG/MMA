#!/bin/bash


for fold in sim_*; do

    # ksimulation=$(echo $simulation | grep -Po '(?<=_)\d+')
    cd $fold
        echo $fold
        mv results_*.h5 ../
        # mv CUPRAD-*.output ../
        # mv pre_processor-*.output ../
    cd ..

done;



echo "done" 
