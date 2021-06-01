#!/bin/bash
#SBATCH -p CMSlab
#SBATCH --job-name=test28
#SBATCH -N 1
#SBATCH -t 11-00:00:00
#SBATCH -n 64
#SBATCH --cpus-per-task=1
#SBATCH --mem=500g
#SBATCH --exclusive
###SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=charles.puelz@gmail.com

##export SLURM_OVERLAP=1
##export MV2_ENABLE_AFFINITY=1
##export MV2_USE_SHARED_MEM=0

MPI=/nas/longleaf/home/mrdavey/heart/programs/mvapich/2.3.3/bin

# $MPI/mpirun ./main3d contraction_tests/input.new_mesh_v150 \
# -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none \
# -velocity_ksp_norm_type none -pressure_ksp_norm_type none \
# -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason -ksp_converged_reason -ksp_monitor_true_residual

#srun --mpi=pmi2 ./main3d contraction_tests/input.new_mesh_v150 -stokes_ksp_monitor -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6 -stokes_ksp_converged_reason -ksp_converged_reason -ksp_monitor_true_residual

srun --mpi=pmi2 ./main3d initial_tests/input.test28 -velocity_ksp_type cg -velocity_pc_type none -velocity_ksp_norm_type none -pressure_ksp_norm_type none -stokes_ksp_rtol 1.0e-6

# /pine/scr/m/r/mrdavey/fung/restart_IB3d 60000 \

