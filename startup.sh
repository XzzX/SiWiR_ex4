#!/bin/bash -l

#PBS -N se_sim
#PBS -l nodes=4:ppn=32
#PBS -l walltime=0:05:00
#PBS -q siwir
#PBS -o $PBS_JOBNAME.out -e $PBS_JOBNAME.err

. /etc/profile.d/modules.sh
module load openmpi/1.6.5-ib
module load gcc/4.8.2

cd ~/SiWiR/siwir_ex3

mpirun -np 128 ./cg 10000 10000 100 -1
