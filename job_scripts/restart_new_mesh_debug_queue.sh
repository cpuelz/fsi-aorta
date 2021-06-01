#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=v17_new_mesh
#SBATCH --ntasks-per-node=32
#SBATCH -N 2
#SBATCH -p debug_queue
#SBATCH -e new_mesh_v17_err
#SBATCH -o new_mesh_v17_out
#SBATCH -t 0-4:0:0
mpirun ./main3d contraction_tests/input.new_mesh_v17 /21dayscratch/scr/c/p/cpuelz/new_mesh_v17/restart_IB3d 110000 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
