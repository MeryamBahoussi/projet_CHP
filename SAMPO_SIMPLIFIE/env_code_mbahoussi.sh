#!/bin/sh

export MACHTYPE
export HOSTARCH=x86_64-intel-linux

#export GCC_VERSION=9.3.0

#export OPENMPI_VERSION=4.0.3

export PETSC_VERSION=3.12.4_1
source ./env_petsc.sh

#export NETCDF_VERSION=4.7.3_2

# pour mutation
export MPP_DIRECTORY=
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
#export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$DYLD_LIBRARY_PATH

