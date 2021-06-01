#!/bin/bash
#SBATCH --exclusive
#SBATCH --job-name=HO_EDPV_test
#SBATCH --ntasks-per-node=16
#SBATCH -N 2
#SBATCH -p skylake
#SBATCH -e HO_EDPV_test_err
#SBATCH -o HO_EDPV_test_out
#SBATCH -t 7-0:0:0
mpirun ./main3d contraction_tests/input.HO_EDPV_test -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason
