#!/bin/sh

#sys="iffslurm"
sys="mac"

make clean
make -f makefile_${sys}            
make clean
make -f makefile_${sys}_py tbfitpy_mod
make clean
make -f makefile_${sys}_mpi_py  tbfitpy_mod
make clean
make -f makefile_${sys}_mpi     
