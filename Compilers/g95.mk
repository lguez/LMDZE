# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For G95 0.91

FC = g95

netcdf_inc_dir = /home/guez_local/include/NetCDF_g95 /home/guez_local/include
netcdf_lib_dir = /home/guez_local/lib /home/guez_local/lib/NetCDF_g95

numer_rec_dir = /home/guez_local/lib/Numer_Rec_Lionel/n
netcdf95_dir = /home/guez_local/lib/NetCDF95/g95
IOIPSL_dir = /home/guez_local/lib/IOIPSL_Lionel/an

# Include flags:
inc_flags = $(addprefix -I, ${libf_dir} ${libf_dir}/dyn3d ${libf_dir}/phylmd ${libf_dir}/filtrez ${netcdf_inc_dir} ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir})

# Other flags which do not affect run time performance:
lang_flags = -ffree-form -pedantic -std=f95 -Wall -Wextra -Wno=136,163,165
# Warning (136): Module variable is never used
# Warning (163): Actual argument does not have an INTENT
# Warning (165): Implicit interface

# Flags which affect run time performance:
perf_flags = -fbounds-check -freal=nan -ftrace=full -g -O0

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir}) -lioipsl -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec
