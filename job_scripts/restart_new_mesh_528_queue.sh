#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=v72_new_mesh
#SBATCH -n 48
#SBATCH -N 2
#SBATCH -p 528_queue
#SBATCH -e new_mesh_v72_err
#SBATCH -o new_mesh_v72_out
#SBATCH -t 3-0:0:0
#SBATCH --mail-user=charles.puelz@gmail.com
#SBATCH --mail-type=ALL

mpirun -n 48 -npernode 24 ./main3d contraction_tests/input.new_mesh_v72 /21dayscratch/scr/c/p/cpuelz/new_mesh_v72/restart_IB3d 006500 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason

