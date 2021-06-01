#!/bin/bash
#SBATCH --job-name=19contraction_test
#SBATCH -n 16
#SBATCH -p debug_queue
#SBATCH -e contraction_test19_err
#SBATCH -o contraction_test19_out
#SBATCH -t 0-4:0:0
mpirun ./main3d contraction_tests/input.contraction19 /21dayscratch/scr/c/p/cpuelz/contraction_test19/restart_IB3d 006000 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
