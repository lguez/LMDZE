# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For G95 0.91

FC = g95

# Include flags:
inc_flags =-I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I/home/guez/include/NetCDF_g95 -I/home/guez/lib/IOIPSL_Lionel/ac -I/home/guez/lib/Numer_Rec_Lionel/a -I/home/guez/lib/NetCDF95_g95

# Other flags which do not affect run time performance:
lang_flags = -ffree-form -pedantic -std=f95 -Wall -Wextra -Wno=136,163,165
# Warning (136): Module variable is never used
# Warning (163): Actual argument does not have an INTENT
# Warning (165): Implicit interface

# Flags which affect run time performance:
perf_flags = -fbounds-check -freal=nan -ftrace=full -g -O0

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDLIBS = -L/home/guez/lib/IOIPSL_Lionel/ac -L/home/guez/lib -L/home/guez/lib/NetCDF_g95 -L/home/guez/lib/Numer_Rec_Lionel/a -L/home/guez/lib/NetCDF95_g95 -lioipsl -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec
