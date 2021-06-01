#!/bin/bash
#SBATCH --job-name=7afterload
#SBATCH -n 16
#SBATCH -p skylake
#SBATCH -e afterload_test7_err
#SBATCH -o afterload_test7_out
#SBATCH -t 7-0:0:0
mpirun ./main3d afterload_tests/input.afterload7 -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
