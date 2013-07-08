# This is a makefile for GNU make.

# This makefile builds LMDZE.

# Suffixes are:
# "f90" for free format, no preprocessing
# "f" for fixed format, no preprocessing

# 1. Source files

src_root = .

VPATH = ${src_root}:${src_root}/dyn3d:${src_root}/dyn3d/Vlsplt:${src_root}/filtrez:${src_root}/phylmd:${src_root}/bibio:${src_root}/phylmd/Mobidic:${src_root}/phylmd/Orography:${src_root}/phylmd/Radlwsw:${src_root}/IOIPSL:${src_root}/IOIPSL/Stringop:${src_root}/dyn3d/Read_reanalyse:${src_root}/phylmd/Thermcell:${src_root}/phylmd/CV3_routines:${src_root}/phylmd/Conflx:${src_root}/phylmd/CV_routines:${src_root}/phylmd/Interface_surf:${src_root}/dyn3d/Dissipation:${src_root}/IOIPSL/Histcom:${src_root}/IOIPSL/Flincom

common_sources := $(shell cat ${src_root}/common_sources)

sources_ce0l := $(sort ${common_sources} $(shell cat ${src_root}/sources_ce0l_only))

sources_gcm := $(sort ${common_sources} $(shell cat ${src_root}/sources_gcm_only))

sources := $(sort ${sources_ce0l} ${sources_gcm})

# 2. Objects and executable files

# Object files:
obj_ce0l := $(addsuffix .o, $(basename ${sources_ce0l}))
obj_gcm := $(addsuffix .o, $(basename ${sources_gcm}))
objects := $(addsuffix .o, $(basename ${sources}))

# Executable files:
execut = ce0l gcm

# 3. Compiler-dependent part

mode = debug
include ../Compilers/${FC}_${mode}.mk

# 4. Rules

SHELL = bash

# Extend known suffixes:
COMPILE.f90 = $(FC) $(F90FLAGS) $(TARGET_ARCH) -c
%.o: %.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) $<

.PHONY: all clean clobber depend

all: ${execut} trace

${execut}:
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

ce0l: ${obj_ce0l}
gcm: ${obj_gcm}

depend ${src_root}/depend.mk:
	makedepf90 -Wmissing -Wconfused -I${VPATH} -nosrc $(addprefix -u, netcdf numer_rec_95 netcdf95 nr_util jumble) ${sources} >${src_root}/depend.mk

TAGS: ${sources}
	ctags -e --language-force=fortran $^

clean:
	rm -f ${objects}
	rm -f ${execut} trace

clobber: clean
	rm -f *.mod ${src_root}/depend.mk ${src_root}/TAGS

trace:
	${FC} ${version_flag} >$@ 2>&1
	echo -e "\nFC = ${FC}\n\nFFLAGS = ${FFLAGS}\n\nLDLIBS = ${LDLIBS}\n\nLDFLAGS = ${LDFLAGS}" >>$@

# Dependencies among object files:
include ${src_root}/depend.mk

# Other rules, optionnally:
-include grep.mk
-include nag_rules.mk
