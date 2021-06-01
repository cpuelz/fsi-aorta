#!/bin/bash
#SBATCH -p skx-normal
#SBATCH --job-name=test3
#SBATCH -N 2
#SBATCH -t 48:00:00
#SBATCH -n 24
#SBATCH -o test3_out
#SBATCH -e test3_err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=charles.puelz@gmail.com

ibrun ./main3d initial_tests/input.test3 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason

