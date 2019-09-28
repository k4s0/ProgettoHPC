#!/bin/bash

rm graph-image/*out*

nsteps=100000
n=256

clear

echo Testing the Serial version...
./earthquake ${nsteps} ${n}
mv serial_out graph-image/

for th in $(seq 1 12);
do
  echo Testing the MPI version with ${th} threads... 
  mpirun -n ${th} mpi-earthquake ${nsteps} ${n} mpi_out_${th}
  cat mpi_out_${th} >> mpi_final
  mv mpi_out_${th} graph-image/
done

rm graph-image/*out*

echo All test are done.

