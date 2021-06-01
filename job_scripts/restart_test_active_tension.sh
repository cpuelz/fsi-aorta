#!/bin/bash
#SBATCH --job-name=test_active_tension
#SBATCH -n 48
#SBATCH -p 528_queue
#SBATCH -e test_active_tension_err
#SBATCH -o test_active_tension_out
#SBATCH -t 3-0:0:0
mpirun --map-by core ./main3d contraction_tests/input.test_active_tension /21dayscratch/scr/c/p/cpuelz/ho_myo/restart_IB3d 19000 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
