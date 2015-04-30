# This is a makefile for GNU make.

# This makefile builds LMDZE.

# 1. Source files

makefile_dir = .

VPATH := ${makefile_dir}/Sources $(addprefix ${makefile_dir}/Sources/, $(shell cat ${makefile_dir}/directories))

src_ce0l := $(shell cat ${makefile_dir}/src_ce0l)
src_gcm := $(shell cat ${makefile_dir}/src_gcm)
src_test_ozonecm := $(shell cat ${makefile_dir}/src_test_ozonecm)
src_test_inter_barxy := $(shell cat ${makefile_dir}/src_test_inter_barxy)
src_test_fxhyp := $(shell cat ${makefile_dir}/src_test_fxhyp)

sources := $(sort ${src_ce0l} ${src_gcm} ${src_test_ozonecm} ${src_test_inter_barxy} ${src_test_fxhyp})

# 2. Objects and executable files

obj_ce0l := $(src_ce0l:.f=.o)
obj_gcm := $(src_gcm:.f=.o)
obj_test_ozonecm := $(src_test_ozonecm:.f=.o)
obj_test_inter_barxy := $(src_test_inter_barxy:.f=.o)
obj_test_fxhyp := $(src_test_fxhyp:.f=.o)
objects := $(sources:.f=.o)
execut = ce0l gcm test_ozonecm test_inter_barxy test_fxhyp

# 3. Compiler-dependent part

mode = debug
include Compiler_options/${FC}_${mode}.mk

# 4. Rules

SHELL = bash
.DELETE_ON_ERROR:
.PHONY: all clean clobber depend
all: ${execut} log

${execut}:
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

ce0l: ${obj_ce0l}
gcm: ${obj_gcm}
test_ozonecm: ${obj_test_ozonecm}
test_inter_barxy: ${obj_test_inter_barxy}
test_fxhyp: ${obj_test_fxhyp}

depend ${makefile_dir}/depend.mk:
	makedepf90 -free -Wmissing -Wconfused $(addprefix -I, ${VPATH}) -nosrc $(addprefix -u , netcdf numer_rec_95 netcdf95 nr_util jumble) ${sources} >${makefile_dir}/depend.mk

TAGS: ${sources}
	ctags -e --language-force=fortran $^

clean:
	rm -f ${execut} ${objects} log

clobber: clean
	rm -f *.mod ${makefile_dir}/depend.mk TAGS

log:
	hostname >$@
	${FC} ${version_flag} >>$@ 2>&1
	echo -e "\nFC = ${FC}\n\nFFLAGS = ${FFLAGS}\n\nLDLIBS = ${LDLIBS}\n\nLDFLAGS = ${LDFLAGS}" >>$@

ifeq ($(findstring $(MAKECMDGOALS), clobber depend),)
include ${makefile_dir}/depend.mk
endif

-include grep.mk
-include nag_rules.mk
