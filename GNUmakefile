# This is a makefile for GNU make.

# This makefile builds LMDZE.

# Suffixes are:
# "f90" for free format, no preprocessing
# "f" for fixed format, no preprocessing

# 1. Source files

src_root = .

VPATH := ${src_root} $(addprefix ${src_root}/, $(shell cat ${src_root}/directories))

common_sources := $(shell cat ${src_root}/common_sources)
src_ce0l_only := $(shell cat ${src_root}/src_ce0l_only)
src_gcm_only := $(shell cat ${src_root}/src_gcm_only)
sources = ${src_ce0l_only} ${src_gcm_only} ${common_sources}

# 2. Objects and executable files

obj_ce0l := $(addsuffix .o, $(sort $(basename ${common_sources} ${src_ce0l_only})))

obj_gcm := $(addsuffix .o, $(sort $(basename ${common_sources} ${src_gcm_only})))

objects := $(addsuffix .o, $(basename ${sources}))
execut = ce0l gcm

# 3. Compiler-dependent part

mode = debug
include Compiler_options/${FC}_${mode}.mk

# 4. Rules

SHELL = bash
COMPILE.f90 = $(FC) $(F90FLAGS) $(TARGET_ARCH) -c

%.o: %.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) $<

.DELETE_ON_ERROR:
.PHONY: all clean clobber depend
all: ${execut} log

${execut}:
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

ce0l: ${obj_ce0l}
gcm: ${obj_gcm}

depend ${src_root}/depend.mk:
	makedepf90 -Wmissing -Wconfused $(addprefix -I, ${VPATH}) -nosrc $(addprefix -u , netcdf numer_rec_95 netcdf95 nr_util jumble) ${sources} >${src_root}/depend.mk

${src_root}/TAGS: ${sources}
	ctags -e --language-force=fortran -f $@ $^

clean:
	rm -f ${execut} ${objects} log

clobber: clean
	rm -f *.mod ${src_root}/depend.mk ${src_root}/TAGS

log:
	hostname >$@
	${FC} ${version_flag} >>$@ 2>&1
	echo -e "\nFC = ${FC}\n\nFFLAGS = ${FFLAGS}\n\nLDLIBS = ${LDLIBS}\n\nLDFLAGS = ${LDFLAGS}" >>$@

ifneq ($(MAKECMDGOALS), clobber)
include ${src_root}/depend.mk
endif

-include grep.mk
-include nag_rules.mk
