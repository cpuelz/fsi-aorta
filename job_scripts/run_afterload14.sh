#!/bin/bash
#SBATCH --job-name=14afterload
#SBATCH --ntasks=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e afterload_test14_err
#SBATCH -o afterload_test14_out
#SBATCH -t 7-0:0:0
mpirun ./main3d afterload_tests/input.afterload14 -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
