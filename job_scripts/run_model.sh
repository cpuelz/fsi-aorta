#!/bin/bash
# make clean && make
mkdir $HOSTNAME
cp *.e $HOSTNAME
cp input.pa_aorta_heart_av $HOSTNAME
cp input.pa_aorta_heart $HOSTNAME
cp input.aorta_heart $HOSTNAME
cp input.heart $HOSTNAME
cd $HOSTNAME
rm nohup.out
nohup mpirun -n 24 --bind-to-core ../main3d ../input.pa_aorta_heart_av -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-6 &
#nohup ../main3d ../input.aorta_and_heart -stokes_ksp_monitor -stokes_ksp_converged_reason -stokes_ksp_rtol 1.0e-8 &
