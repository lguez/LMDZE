# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

FC = gfortran

netcdf_inc_dir = /user/guez_local/include
netcdf_lib_dir = /user/guez_local/lib

numer_rec_dir = /user/guez_local/lib/Numer_Rec_Lionel_b
nr_util_dir = /user/guez_local/lib/NR_util_j
netcdf95_dir = /user/guez_local/lib/NetCDF95_gfortran

# Include flags:
inc_flags = $(addprefix -I, ${libf_dir} ${libf_dir}/phylmd ${netcdf_inc_dir} ${numer_rec_dir} ${netcdf95_dir} ${nr_util_dir})

# Other flags which do not affect run time performance:
lang_flags = -ffree-form -frange-check -std=f95 -pedantic-errors -Wall -Wconversion -Wimplicit-interface -Wunderflow -Wextra -Wunreachable-code

# Flags which affect run time performance:
perf_flags = -fbacktrace -ffpe-trap=invalid,zero,overflow -fbounds-check -g3 -O0 -fstack-protector-all

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_dir} ${netcdf95_dir} ${nr_util_dir}) -lnetcdf95 -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lnumer_rec -lnr_util
