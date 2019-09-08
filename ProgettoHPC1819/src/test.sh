#!/bin/bash

nsteps=100000
n=256

echo > graph-image/serial
./earthquake ${nsteps} ${n} > graph-image/serial

for th in $(seq 1 12);
do
  echo > graph-image/omp-${th} 
  OMP_NUM_THREADS=${th} ./omp-earthquake ${nsteps} ${n} > graph-image/omp-${th}
done

make clean
