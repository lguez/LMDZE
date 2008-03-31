# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For G95 0.91

FC = g95

netcdf_inc_dir = /home/guez/include/NetCDF_g95
netcdf_lib_dir = /home/guez/lib /home/guez/lib/NetCDF_g95

numer_rec_dir = /home/guez/lib/Numer_Rec_Lionel/a
netcdf95_dir = /home/guez/lib/NetCDF95/g95
IOIPSL_dir = /home/guez/lib/IOIPSL_Lionel/ac

# Include flags:
inc_flags =-I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I${netcdf_inc_dir} -I${numer_rec_dir} -I${netcdf95_dir} -I${IOIPSL_dir}

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
