#!/bin/bash
#SBATCH --job-name=10contraction_test
#SBATCH -n 16
#SBATCH -p skylake
#SBATCH -e contraction_test10_err
#SBATCH -o contraction_test10_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.contraction10 -stokes_ksp_monitor -stokes_ksp_converged_reason
