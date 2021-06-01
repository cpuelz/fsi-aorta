#!/bin/bash
#SBATCH --job-name=pap_contraction_test
#SBATCH --ntasks=8
#SBATCH -N 1
#SBATCH -p debug_queue
#SBATCH -e pap_contraction_test_err
#SBATCH -o pap_contraction_test_out
#SBATCH -t 0-4:0:0
mpirun ./main3d contraction_tests/input.pap_contraction_test -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
