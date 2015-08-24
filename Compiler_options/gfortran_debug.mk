# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /usr/include
netcdf_lib_dir = 

numer_rec_95_dir = ${HOME}/Desktop/lib/Numer_Rec_95_debug
nr_util_dir = ${HOME}/Desktop/lib/NR_util_debug
netcdf95_dir = ${HOME}/Desktop/lib/NetCDF95_debug
jumble_dir = ${HOME}/Desktop/lib/Jumble_debug

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Fortran language options:
FFLAGS += -ffree-form -std=f95

# Error and warning options:
FFLAGS += -fmax-errors=1 -pedantic-errors -Wall -Wcharacter-truncation -Wunderflow -Wunreachable-code -Wno-conversion

# Debugging options:
FFLAGS += -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -g

# Code generation options:
FFLAGS += -fcheck=bounds -fcheck=do -fcheck=mem -fcheck=pointer -fcheck=recursion -finit-real=SNAN

# Optimization options:
FFLAGS += -O0

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}) -ljumble -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec_95 -lnr_util

version_flag = --version