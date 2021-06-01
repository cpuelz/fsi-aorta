#!/bin/bash
#SBATCH --job-name=test5
#SBATCH --ntasks=32
#SBATCH --mem-per-cpu=1GB
#SBATCH --constraint=opath
#SBATCH --time=24:00:00
#SBATCH --mail-user=cp16@rice.edu
#SBATCH --mail-type=ALL

srun main3d initial_tests/input.test5 -ksp_converged_reason -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
