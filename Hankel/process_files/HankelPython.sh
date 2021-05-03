#!/bin/bash

module purge
module load python/3.7.2

cp $HOME/git/Hankel/Hankel.py .
cp $HOME/git/Hankel/Hankel1.py .
cp $HOME/git/Hankel/Hankelserial.py .

cd res1
  rm Spectrum.dat
  rm Spectrumreal.dat
  rm Spectrumimag.dat
  rm omegagrid_anal.dat
  rm rgrid_anal.dat
cd ..

#python Hankel1.py
python Hankel.py
#python Hankelserial.py
