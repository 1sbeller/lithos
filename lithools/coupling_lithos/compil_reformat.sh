#!/bin/bash

F90='mpif90 -O3 -ipo -xHost -assume byterecl'

$F90 -c ../../common/precision_mod.f90
$F90 -c ../../common/constants_mod.f90
$F90 -c gll_library.f90
$F90 -c interpolation_parameters_mod.f90
$F90 -c mpi_int_mod.f90
$F90 -c inputs_outputs_interp_mod.f90
$F90 -c interpolate_3D_wavefield.f90
$F90 *.o -o interpolate_3D_wavefield.x
 

