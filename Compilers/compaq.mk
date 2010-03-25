# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For the MIPSpro 7 Fortran 90 compiler

FC = f90

# Include flags:
inc_flags = -I${libf_dir} -I${libf_dir}/phylmd -I/usr/local/pub/include -I${workdir}/IOIPSL_k

# Other flags which do not affect run time performance:
lang_flags = -ansi -fullwarn

# Flags which affect run time performance:
perf_flags = -check_bounds -g2 -O0 -DEBUG:div_check=3:subscript_check=ON:verbose_runtime=ON -DEBUG:trap_uninitialized=ON

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}
LDFLAGS =

LDLIBS= -L${workdir}/IOIPSL_k -lioipsl -L/usr/local/pub/lib64 -lnetcdf
