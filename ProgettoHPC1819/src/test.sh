#!/bin/bash

rm graph-image/*_out
rm *_out

touch serial_out
touch omp_out

nsteps=100000
n=256

clear

echo Testing the serial version...
./earthquake ${nsteps} ${n} > serial_out
mv serial_out graph-image/

for th in $(seq 1 12);
do
  echo Testing the OpenMP version with ${th} threads... 
  OMP_NUM_THREADS=${th} ./omp-earthquake ${nsteps} ${n} >> omp_out 
  mv omp_out graph-image/ 
done

echo All test are done.

