# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /usr/include
netcdf_lib_dir = 

numer_rec_95_dir = /user/guez_local/lib/Numer_Rec_95_gfortran_debug
nr_util_dir = /user/guez_local/lib/NR_util_gfortran_debug
netcdf95_dir = /user/guez_local/lib/NetCDF95_gfortran_debug
jumble_dir = /user/guez_local/lib/Jumble_debug

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Fortran language options:
FFLAGS += -std=f95

# Error and warning options:
FFLAGS += -fmax-errors=1 -pedantic-errors -Wall -Wcharacter-truncation -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code

# Debugging options:
FFLAGS += -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -g

# Code generation options:
FFLAGS += -fcheck=all -finit-real=SNAN

# Optimization options:
FFLAGS += -O0

F90FLAGS = ${FFLAGS}

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}) -ljumble -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec_95 -lnr_util

version_flag = --version
