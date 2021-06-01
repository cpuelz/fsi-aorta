#!/bin/bash
#SBATCH --job-name=prestrain
#SBATCH --ntasks=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e prestrain_err
#SBATCH -o prestrain_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.prestrain_test -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
