#!/bin/bash

rm graph-image/*out*

nsteps=100000
n=256

clear

echo Testing the serial version...
./earthquake ${nsteps} ${n}
mv serial_out graph-image/
echo
for th in $(seq 1 12);
do
  echo Testing the OpenMP version with ${th} threads... 
  OMP_NUM_THREADS=${th} ./omp-earthquake ${nsteps} ${n} omp_out_${th} 
  cat omp_out_${th} >> omp_final
  mv omp_out_${th} graph-image/
done
echo
for th in $(seq 1 12);
do
  echo Testing the MPI version with ${th} threads... 
  mpirun -n ${th} mpi-earthquake ${nsteps} ${n} mpi_out_${th}
  cat mpi_out_${th} >> mpi_final
  mv mpi_out_${th} graph-image/
done

rm graph-image/*out*

echo All test are done.

