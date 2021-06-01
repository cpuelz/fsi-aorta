#!/bin/bash
#SBATCH --job-name=6contraction_test
#SBATCH -n 24
#SBATCH -p skylake
#SBATCH -e contraction_test6_err
#SBATCH -o contraction_test6_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.contraction6 -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
