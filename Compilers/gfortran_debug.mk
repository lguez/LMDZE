# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /usr/include
netcdf_lib_dir = 

numer_rec_95_dir = /user/guez_local/lib/Numer_Rec_95_gfortran_debug
nr_util_dir = /user/guez_local/lib/NR_util_gfortran_debug
netcdf95_dir = /user/guez_local/lib/NetCDF95_gfortran_debug
jumble_dir = /user/guez_local/lib/Jumble_debug

# Include flags:
inc_flags = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Other flags which do not affect run time performance:
lang_flags = -ffree-form -frange-check -std=f95 -pedantic-errors -Wall -Wunderflow -Wextra

# Flags which affect run time performance:
perf_flags = -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}) -ljumble -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec_95 -lnr_util

version_flag = --version
