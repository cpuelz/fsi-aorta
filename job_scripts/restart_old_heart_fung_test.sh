#!/bin/bash
#SBATCH --job-name=fung_test
#SBATCH --ntasks=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e old_heart_fung_test_err
#SBATCH -o old_heart_fung_test_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.old_heart_fung_test /21dayscratch/scr/c/p/cpuelz/old_heart_fung_test/restart_IB3d 084000 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
