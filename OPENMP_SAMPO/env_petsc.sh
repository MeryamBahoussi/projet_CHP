#!/bin/sh

export PETSC_DIR=/usr/local/opt/petsc
export PETSC_ARCH=
export PETSC_LIB=${PETSC_DIR}/${PETSC_ARCH}/lib
export PETSC_BIN=${PETSC_DIR}/${PETSC_ARCH}/bin
export PETSC_INC=${PETSC_DIR}/${PETSC_ARCH}/include
export PETSC_LDFLAGS=-L${PETSC_DIR}/${PETSC_ARCH}/lib
