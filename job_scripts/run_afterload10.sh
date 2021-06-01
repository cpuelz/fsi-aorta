#!/bin/bash
#SBATCH --job-name=10afterload
#SBATCH -n 16
#SBATCH -p skylake
#SBATCH -e afterload_test10_err
#SBATCH -o afterload_test10_out
#SBATCH -t 7-0:0:0
mpirun ./main3d afterload_tests/input.afterload10 -build_twosided allreduce -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
