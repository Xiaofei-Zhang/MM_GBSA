#!/bin/bash

#SBATCH -p Long
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --mail-type=ALL

source /usr/share/Modules/init/bash
module load amber/openmpi/12p12

mpiexec -np 16 python mmgbsa_MPI.py lig rcptr 1.conf
