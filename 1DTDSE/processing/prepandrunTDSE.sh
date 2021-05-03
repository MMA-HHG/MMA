#!/bin/bash

# compile
cd programs
  programspath=`pwd`
  ./compile1DTDSE.sh &> compilationlog.log
  echo "compiled"
  mv compilationlog.log ../logfiles/1DTDSEcompilationlog.log
cd ..

#distribute & submit on cluster

echo $programspath


   cd sims
#   read PRINTFIELDS < FLAGS;
#   if (( $PRINTFIELDS == 1 )); then

     for fold1 in z_*; do

#       echo $fold1
       cd $fold1
       for fold2 in r_*; do
#         echo $fold2
         cd $fold2
           cp $programspath/TDSE1D.out .
           cp $programspath/submit_1DTDSE.slurm .
           jobname="${fold2}_${fold1}"
           echo $jobname
#           sbatch -J TDSE_$jobname submit_1DTDSE.slurm
         cd ..
       done

       cd ..
#       fold3=${fold2:6};
#       fold3=${fold3%.dat};
#       echo $fold3
#       cd ../TDSEs
#       mkdir $fold3
#       cd $fold3
#       mkdir results
#       mkdir inputs
#       # copy neccessary files
#       cp ../../fields/$fold2 inputs/Field.dat
#       cp ../../fields/tgrid.dat inputs/tgrid.dat
#      
#       cp ../../../programs/TDSE1D.out TDSE1D.out
#       cp ../../../programs/param.txt param.txt
#       cp ../../../programs/submit_1DTDSE.slurm submit_1DTDSE.slurm

##       ./run1DTDSE
#       jobname="${fold}_${fold3}"
#       echo $jobname
##       qsub -N 1DTDSE_$jobname submit_1DTDSE
##       sbatch -J TDSE_z_{$kdum}_r_${k2} submit_1DTDSE.slurm
##       sbatch -J TDSE_$jobname submit_1DTDSE.slurm



       
     done

#     cp $fold/figures/propagation.jpg figures/propagation_$fold.jpg
#     cp $fold/figures/propagation.pdf figures/propagation_$fold.pdf
     cd ..
#   fi
   cd ..




#cd programs/TDSE
##remove old programs
#rm *.out

#gfortran -c recalcul_params.f90
#gfortran -o calculator.out *.o

#rm *.o

#gfortran -c convert_inps.f90
#gfortran -o convert.out *.o

#rm *.o

#echo "compiled"
#cd ..
#cp programs/*.out .
#cp programs/script2.sh .

##mkdir intensity_study

## clearing folder
#rm -r sim*
#rm params_check.txt
#rm ListOfSimulations.txt
#rm -r logfiles

#mkdir logfiles

## inputs (I0 [W/cm2], press [bar], w0 [SI], focus [SI], tau [fs], n2atm [cm2/W], lambda [SI], THEORY, medium_length [SI]), put a negative number for no lense
## [1-Euler method / 2- Semi-analytical / 3- PPT / 4-ADK-molecular / 5- Keldysh for crystal only / 6- KDP Guillaume / 7- air Pavel / 8-CR]

##run=(1 2 3 4 5 6 7);
##run=(1 2 3 4 5 6); I0=${run[0]}; press=${run[1]}; w0=${run[2]}; tau=${run[3]}; n2atm=${run[4]}; lambda=${run[5]}; export I0 press w0 tau n2atm lambda j; ./script2.sh; let j++;

#j=1
#while read LINE; do 
#	read -a run <<< $LINE;
#	echo "${run[*]}" >> "params_check.txt";
#	I0=${run[0]}; press=${run[1]}; w0=${run[2]}; focus=${run[3]}; tau=${run[4]}; n2atm=${run[5]}; lambda=${run[6]}; THEORY=${run[7]}; medium_length=${run[8]}; PRINTFIELDS=${run[9]};
#	export I0 press w0 focus tau n2atm lambda j THEORY medium_length PRINTFIELDS; ./script2.sh; let j++; 
#done < params.txt


##run=(3 4); let j++; n2atm=${run[0]}; n2atm=n2atm/2.0 | bc; export n2atm j; ./script2.sh;

#echo "test4"


#for fold in sim*; do echo $fold; echo $fold >> "ListOfSimulations.txt"; done;

## clearing folder

#mv ListOfSimulations.txt logfiles/
#mv params_check.txt logfiles/
#mv *.dat logfiles/
#rm *.out
#rm script2.sh

##if (( $(echo "$focus < 0" | bc -l) )); then
##  echo "no lense"
##fi



echo "done" 
