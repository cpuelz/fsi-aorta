#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=stronger_first_order_ho_myo
#SBATCH --ntasks-per-node=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e ho_myo_first_order_stronger_err
#SBATCH -o ho_myo_first_order_stronger_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.ho_myo_first_order_stronger_contraction -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
