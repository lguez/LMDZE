# This is the makefile for "test_ozonecm". For GNU Make 3.80.

# Extend known suffixes:
%.o: %.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
COMPILE.f90 = $(FC) $(F90FLAGS) -c

# Directories:
srcdir = .

VPATH = ${srcdir}:${srcdir}/../phylmd:${srcdir}/../dyn3d

# Compiler-dependent part:
FC = g95

netcdf_inc_dir = /home/guez_local/include/NetCDF_g95 /home/guez_local/include
netcdf_lib_dir = /home/guez_local/lib/NetCDF_g95 /home/guez_local/lib

numer_rec_dir = /home/guez_local/lib/Numer_Rec_Lionel/r
netcdf95_dir = /home/guez_local/lib/NetCDF95/g95

# Include flags:
inc_flags = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_dir} ${netcdf95_dir})

# Other flags which do not affect run time performance:
lang_flags = -ffree-form -pedantic -std=f95 -Wall -Wextra

# Flags which affect run time performance:
perf_flags = -fbounds-check -freal=nan -ftrace=full -g -O0

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_dir} ${netcdf95_dir}) -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec

sources = test_ozonecm.f90 dimens_m_test.f90 comvert.f90 ozonecm.f90 phyetat0_test.f90 dimphy_test.f90

objects = $(sources:.f90=.o)

# Executable file:
execut = test_ozonecm

${execut}: ${objects}
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

TAGS: ${sources}
	etags -l fortran $^

.PHONY: clean clobber

clean:
	rm -f ${objects}
	rm -f ${execut}

clobber: clean
	rm -f *.mod ${srcdir}/depend.mk TAGS

.PHONY: depend
depend ${srcdir}/depend.mk:
	makedepf90 -Wmissing -Wconfused -I ${VPATH} -u netcdf -u numer_rec -u netcdf95 -nosrc ${sources} >${srcdir}/depend.mk

# Dependencies between object files and include files:
include ${srcdir}/depend.mk
