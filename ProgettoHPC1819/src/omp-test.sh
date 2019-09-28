#!/bin/bash

rm graph-image/*out*

nsteps=100000
n=256

clear

for th in $(seq 1 12);
do
  echo Testing the OpenMP version with ${th} threads... 
  OMP_NUM_THREADS=${th} ./omp-earthquake ${nsteps} ${n} omp_out_${th} 
  cat omp_out_${th} >> omp_final
  mv omp_out_${th} graph-image/
done

rm graph-image/*out*

echo All test are done.

