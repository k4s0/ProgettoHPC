#!/bin/bash

rm graph-image/*out*
rm serial_out
rm omp_final

touch serial_out
touch omp_out

nsteps=100000
n=256

clear

echo Testing the serial version...
./earthquake ${nsteps} ${n}
mv serial_out graph-image/

for th in $(seq 1 12);
do
  echo Testing the OpenMP version with ${th} threads... 
  OMP_NUM_THREADS=${th} ./omp-earthquake ${nsteps} ${n} omp_out_${th} 
  cat omp_out_${th} >> omp_final
  mv omp_out_${th} graph-image/ 
done

echo All test are done.

