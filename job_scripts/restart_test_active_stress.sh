#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=test_active_stress
#SBATCH --ntasks-per-node=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e test_active_stress_err
#SBATCH -o test_active_stress_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.test_active_stress /21dayscratch/scr/c/p/cpuelz/test_active_stress/restart_IB3d 69000 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
