#!/bin/bash
#SBATCH --job-name=whole_heart_simulation
#SBATCH -n 64
#SBATCH -p 528_queue
#SBATCH -e err
#SBATCH -o out
#SBATCH -t 3-0:0:0
mpirun ./main3d input.test15 -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
