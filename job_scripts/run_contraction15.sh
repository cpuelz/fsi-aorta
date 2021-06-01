#!/bin/bash
#SBATCH --job-name=15contraction_test
#SBATCH -n 24
#SBATCH -p skylake
#SBATCH -e contraction_test15_err
#SBATCH -o contraction_test15_out
#SBATCH -t 7-0:0:0
mpirun --map-by core ./main3d contraction_tests/input.contraction15 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
