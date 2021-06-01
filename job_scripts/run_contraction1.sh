#!/bin/bash
#SBATCH --job-name=1contraction_test
#SBATCH -n 16
#SBATCH -p knl
#SBATCH -e contraction_test1_err
#SBATCH -o contraction_test1_out
#SBATCH -t 2-0:0:0
mpirun ./main3d contraction_tests/input.Guccione -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6
