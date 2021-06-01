#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=1st_av_test
#SBATCH --ntasks-per-node=32
#SBATCH -N 2
#SBATCH -p 528_queue
#SBATCH -e 1st_av_test_err
#SBATCH -o 1st_av_test_out
#SBATCH -t 3-0:0:0
mpirun ./main3d contraction_tests/input.av_test_first -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
