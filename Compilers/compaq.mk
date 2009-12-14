# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For the MIPSpro 7 Fortran 90 compiler on Rhodes.

FC = f90

FFLAGS = -I${bypr_dir} -Igrid -Idyn3d -Iphylmd -I/usr/local/pub/include -I${workdir}/IOIPSL_k

F90FLAGS = -ansi -check_bounds -fullwarn -g2 -O0 -DEBUG:div_check=3:subscript_check=ON:verbose_runtime=ON -DEBUG:trap_uninitialized=ON ${FFLAGS}

CPPFLAGS = $(addprefix -D, ${macros})
LDFLAGS =

LDLIBS=-L${workdir}/IOIPSL_k -lioipsl -L/usr/local/pub/lib64 -lnetcdf
