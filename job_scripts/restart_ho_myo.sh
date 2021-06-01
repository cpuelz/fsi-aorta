#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=ho_myo
#SBATCH --ntasks-per-node=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e ho_myo_err
#SBATCH -o ho_myo_out
#SBATCH -t 7-0:0:0
mpirun --map-by core ./main3d contraction_tests/input.ho_myo /21dayscratch/scr/c/p/cpuelz/ho_myo/restart_IB3d 110000 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
