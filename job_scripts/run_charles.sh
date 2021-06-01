#!/bin/bash
#SBATCH --job-name=charles
#SBATCH --ntasks-per-node=18
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e charles_err
#SBATCH -o charles_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.charles -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
