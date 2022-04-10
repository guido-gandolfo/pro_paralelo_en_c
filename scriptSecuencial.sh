#!/bin/bash
#SBATCH -N 1

#SBATCH -o /nethome/sdyp16/Salida/output.txt
#SBATCH -e /nethome/sdyp16/Errores/errores.txt
mpirun  modSecuencial $1
