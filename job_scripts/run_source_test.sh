#!/bin/bash
#SBATCH --job-name=source_test
#SBATCH -n 16
#SBATCH -p skylake
#SBATCH -e source_err
#SBATCH -o source_out
#SBATCH -t 7-0:0:0
mpirun ./main3d misc_tests/input.sourcetest -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
