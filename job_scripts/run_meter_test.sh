#!/bin/bash
#SBATCH --job-name=metertest
#SBATCH -n 48
#SBATCH -p 528_queue
#SBATCH -e err
#SBATCH -o out
#SBATCH -t 3-0:0:0
mpirun ./main3d input.metertest -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
