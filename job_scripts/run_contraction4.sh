#!/bin/bash
#SBATCH --job-name=4contraction_test
#SBATCH -n 24
#SBATCH -p skylake
#SBATCH -e contraction_test4_err
#SBATCH -o contraction_test4_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.contraction4 -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
