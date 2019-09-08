#!/bin/bash

rm graph-image/omp
rm graph-image/serial

nsteps=100000
n=256

clear

echo Testing the serial version...
echo > graph-image/serial
./earthquake ${nsteps} ${n} > graph-image/serial

for th in $(seq 1 12);
do
  echo Testing the OpenMP versione with ${th} threads... 
  OMP_NUM_THREADS=${th} ./omp-earthquake ${nsteps} ${n} >> graph-image/omp
done

echo All test are done.

make clean
