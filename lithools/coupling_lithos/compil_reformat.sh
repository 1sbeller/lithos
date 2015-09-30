#!/bin/bash

clear
module load intel/15.0.0.090
#module load intelmpi/5.0.1.035
#module load bullmpi/1.2.8.1
module load bullxmpi/1.2.8.3
#source obj/tosource_on.sh
#source /softs/env_default.sh 
 
./clean.sh

F90='mpif90 -O3 -xHost -assume byterecl -traceback' #  -check all -g -warn all'

$F90 -c ../../common/precision_mod.f90
$F90 -c ../../common/constants_mod.f90
$F90 -c interpolation_parameters_mod.f90
$F90 -c mpi_int_mod.f90
$F90 -c interp_process_mod.f90
$F90 -c inputs_outputs_interp_mod.f90
$F90 -c interpolate_3D_wavefield.f90
$F90 *.o -o interpolate_3D_wavefield.x
 
./clean.sh

$F90 -c ../../common/precision_mod.f90
$F90 -c ../../common/constants_mod.f90
$F90 -c gll_library.f90
$F90 -c global_parameters_mod.f90
$F90 -c rotation_matrix_mod.f90
$F90 -c mpi_mod.f90
$F90 -c inputs_outputs_mod.f90
$F90 -c interp_mesh_mod.f90
$F90 -c post_process_wavefields_mod.f90
$F90 -c reconstruct_3D_wavefield.f90
$F90 *.o -o reconstruct_3D_wavefield.x

