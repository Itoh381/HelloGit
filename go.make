#!/bin/bash
export FOPT="-O3"

#FCC=ifort
FCC=gfortran
#/bin/rm -f *.o *.mod
${FCC} -c ${FOPT} precision.f90
${FCC} -c ${FOPT} d3q15model.f90
${FCC} -c ${FOPT} param.f90
${FCC} -c ${FOPT} cavity3d_new_param.f90
${FCC} -c ${FOPT} cavity3d_new.f90
${FCC}    ${FOPT} -o ./lbm *.o

#set FOPT = "-C -traceback"




