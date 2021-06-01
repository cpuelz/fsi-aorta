#!/bin/bash
# Number of nodes
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH -p debug_queue
#SBATCH --time=4:00:00
#SBATCH --job-name=nodes-8-heart-model
#SBATCH --output=nodes-8-heart-model.log


# Load the compiler and MPI library I compiled the program with
module load openmpi_gcc/4.8.5

# Total number of MPI tasks will be calculated by slurm based on either the defaults or command line parameters.

# Need to export the paths to all the HPC Toolkit binaries

# HPCSTRUCT is the tool that you use to do the static analysis of the compiled binary
# This generates information about all the source code that was used to create the application
# and is merged with the measurement data to generate the database that "hpcviewer" and "hpctraceviewer" can open
export HPCSTRUCT=/nas/longleaf/home/deleeke/dogwood/sfw/hpctoolkit/INSTALL/bin/hpcstruct

# HPCRUN is the wrapper that you run the program through to generate the measurements of events at runtime
export HPCRUN=/nas/longleaf/home/deleeke/dogwood/sfw/hpctoolkit/INSTALL/bin/hpcrun

# HPCPROF is used after the measurments and static analysis have been created. It mereges the two datasets to associate the measurment data with the code that was executed.
export HPCPROF=/nas/longleaf/home/deleeke/dogwood/sfw/hpctoolkit/INSTALL/bin/hpcprof-mpi
echo "PID IS $$"

# I export the name of the compiled program to a an environment variable because it is used in the names of several files and directories created by the hpctoolkit tools. This way you can just update the reference once.
export PROGNAME=main3d
export PROG=$PWD/$PROGNAME
export INPUT=$PWD/input.test13

export ROOTDIR=$PWD
$HPCSTRUCT $PROG
export STRUCT=$ROOTDIR/$PROGNAME.hpcstruct

# create a directory to dump all the measurement data and job information into
mkdir experiment
cd experiment

ln -s ${PROG} $PROGNAME
cp $INPUT input3d
ln -s $STRUCT $PROGNAME.hpcstruct
cp $ROOTDIR/*.e .

# Log what the job ID was and what nodes we were using
# So we can reference infromation available from scontrol and sinfo later
touch README
echo "SLURM_JOBID=${SLURM_JOBID}" >> README
echo "SLURM_JOB_NAME=${SLURM_JOB_NAME}" >> README
echo "SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR}" >> README
echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST}" >> README
echo "SLURM_JOB_NUM_NODES=${SLURM_JOB_NUM_NODES}" >> README
echo "SLURM_NTASKS=${SLURM_NTASKS}" >> README

time mpirun $HPCRUN --event REALTIME@6000 ./$PROGNAME input3d 2>&1 | tee -a output.txt

# This is the call that merges the static analysis data and the runtime measurments
$HPCPROF -S $PROGNAME.hpcstruct -I./'*' hpctoolkit-$PROGNAME-measurements* -M stats
