# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /usr/include
netcdf_lib_dir = 

numer_rec_95_dir = ${HOME}/Desktop/lib/Numer_Rec_95_debug
nr_util_dir = ${HOME}/Desktop/lib/NR_util_debug
netcdf95_dir = ${HOME}/Desktop/lib/NetCDF95_debug
jumble_dir = ${HOME}/Desktop/lib/Jumble_debug

lib_dir_list = ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Fortran language options:
FFLAGS += -ffree-form -std=f95

# Error and warning options:
FFLAGS += -fmax-errors=1 -pedantic-errors -Wall -Wcharacter-truncation -Wunderflow -Wunreachable-code -Wno-conversion

# Debugging options:
FFLAGS += -ffpe-trap=invalid,zero,overflow -fbacktrace -fdump-core -g

# Code generation options:
FFLAGS += -fcheck=bounds -fcheck=do -fcheck=mem -fcheck=pointer -fcheck=recursion
##FFLAGS += -finit-real=nan
FFLAGS += -finit-real=SNAN

# Optimization options:
FFLAGS += -O0

# Hardware model options:
FFLAGS += -mcmodel=medium

comma = ,

LDLIBS = $(addprefix -L, ${lib_dir_list}) -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec_95 -ljumble -lnr_util $(addprefix -Wl${comma}-rpath${comma}, ${lib_dir_list})

version_flag = --version
