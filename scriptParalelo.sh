#!/bin/bash
#SBATCH -N 2
#SBATCH --exclusive
#SBATCH --tasks-per-node=1
#SBATCH -o /nethome/sdyp16/Salida/output.txt
#SBATCH -e /nethome/sdyp16/Errores/errores.txt
export OMP_NUM_THREADS=8
mpirun --bind-to none modParalel $1 $2
