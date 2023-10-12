#!/bin/sh

#PBS -q sxs
#PBS --venode  16
#PBS -l elapstim_req=100:00:00
cd $PBS_O_WORKDIR

export VE_PROGINF=DETAIL

source /opt/nec/ve/mpi/3.4.0/bin/necmpivars.sh

mpirun -np 128 ./LetsDoThis 2> proginf.txt &

