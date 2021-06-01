#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=refined_new_myo_parms
#SBATCH --ntasks-per-node=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e refined_new_myo_parms_err
#SBATCH -o refined_new_myo_parms_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.new_myo_parms_refined -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
