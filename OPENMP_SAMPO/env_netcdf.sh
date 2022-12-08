#!/bin/sh

export NETCDF_DIR=/usr/local/opt/netcdf
export NETCDF_ARCH=
export NETCDF_LIB=${NETCDF_DIR}/${NETCDF_ARCH}/lib
export NETCDF_BIN=${NETCDF_DIR}/${NETCDF_ARCH}/bin
export NETCDF_INC=${NETCDF_DIR}/${NETCDF_ARCH}/include
export NETCDF_LDFLAGS=-L${NETCDF_DIR}/${NETCDF_ARCH}/lib
