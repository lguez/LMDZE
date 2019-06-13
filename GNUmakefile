# This is a makefile for GNU make.

# This makefile builds LMDZE.

# 1. Source files

makefile_dir = .

src_ce0l := $(shell cat ${makefile_dir}/src_ce0l)
src_gcm := $(shell cat ${makefile_dir}/src_gcm)
src_test_ozonecm := $(shell cat ${makefile_dir}/src_test_ozonecm)
src_test_inter_barxy := $(shell cat ${makefile_dir}/src_test_inter_barxy)
src_test_fxhyp := $(shell cat ${makefile_dir}/src_test_fxhyp)
src_test_inifilr := $(shell cat ${makefile_dir}/src_test_inifilr)
src_test_orbite = test_orbite.f90 orbite.f90 YOMCST.f90 unit_nml_m.f90

sources := $(sort ${src_ce0l} ${src_gcm} ${src_test_ozonecm} ${src_test_inter_barxy} ${src_test_fxhyp} ${src_test_inifilr} ${src_test_orbite})

netcdf95_dir = ${HOME}/build/Libraries_debug/NetCDF95
nr_util_dir = ${HOME}/build/Libraries_debug/NR_util
numer_rec_95_dir = ${HOME}/build/Libraries_debug/Numer_Rec_95
jumble_dir = ${HOME}/build/Libraries_debug/Jumble

VPATH := ${makefile_dir}

# 2. Compiler-dependent part

mode = debug

comma:= ,
empty:=
space:= $(empty) $(empty)
CPPFLAGS = $(addprefix -D, $(subst ${comma},${space},${cpp_macros}))

inc_dir_list = ${HOME}/build/Libraries_debug/modules /usr/include

lib_dir_list = ${contour_531_dir} ${fortrangis_dir} ${geometry_dir} ${gpc_f_dir} ${jumble_dir} ${netcdf95_dir} ${nr_util_dir} ${numer_rec_95_dir} ${shapelib_03_dir} ${water_dir} ${HOME}/.local/lib

# Include flags:
FFLAGS := $(addprefix -I, ${inc_dir_list})

# Fortran language options:
FFLAGS += -std=f2003

# Error and warning options:
FFLAGS += -fmax-errors=1 -pedantic -Wall -Wcharacter-truncation -Wunused-parameter -Wno-conversion -Wno-integer-division

# Debugging options:
FFLAGS += -fbacktrace -g -ffpe-trap=invalid,zero,overflow

# Code generation options:
FFLAGS += -fcheck=bounds,do,mem,pointer,recursion -fpic -finit-real=nan

# Optimization options:
FFLAGS += -O0

LDLIBS = $(addprefix -L, ${lib_dir_list}) $(addprefix -l, ${lib_list})
version_flag = --version

# 3. General rules

SHELL = bash

%: %.o
	@echo "Linking $@..."
	$(LINK.o) $^ $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.f90
	@echo "Building $@..."
	$(COMPILE.f) $(OUTPUT_OPTION) $<

LINK.o = $(FC) $(LDFLAGS) $(TARGET_ARCH)
.DELETE_ON_ERROR:
.PHONY: all clean clobber depend
all: log

TAGS: ${sources}
	ctags -e --language-force=fortran $^

log:
	hostname >$@
	${FC} ${version_flag} >>$@ 2>&1
	@echo -e "\nFC = ${FC}\n\nCPPFLAGS = ${CPPFLAGS}\n\nFFLAGS = ${FFLAGS}\n\nLDLIBS = ${LDLIBS}\n\nLDFLAGS = ${LDFLAGS}\n\nldflags_lib_dyn = ${ldflags_lib_dyn}" >>$@

clobber: clean
	rm -f *.mod ${makefile_dir}/depend.mk TAGS

VPATH += $(addprefix ${makefile_dir}/, $(shell cat ${makefile_dir}/directories))
cpp_macros = CPP_IIM=16,CPP_JJM=12,CPP_LLM=11
lib_list = numer_rec_95 jumble nr_util netcdf95 netcdff

# 4. Objects and executable files

obj_ce0l := $(addsuffix .o, $(basename ${src_ce0l}))
obj_gcm := $(addsuffix .o, $(basename ${src_gcm}))
obj_test_ozonecm := $(addsuffix .o, $(basename ${src_test_ozonecm}))
obj_test_inter_barxy := $(addsuffix .o, $(basename ${src_test_inter_barxy}))
obj_test_fxhyp := $(addsuffix .o, $(basename ${src_test_fxhyp}))
obj_test_inifilr := $(addsuffix .o, $(basename ${src_test_inifilr}))
obj_test_orbite := $(addsuffix .o, $(basename ${src_test_orbite}))
objects := $(addsuffix .o, $(basename ${sources}))
execut = ce0l gcm test_ozonecm test_inter_barxy test_fxhyp test_inifilr test_orbite

# 5. Specific rules

all: ${execut} log
ce0l: ${obj_ce0l}
gcm: ${obj_gcm}
test_ozonecm: ${obj_test_ozonecm}
test_inter_barxy: ${obj_test_inter_barxy}
test_fxhyp: ${obj_test_fxhyp}
test_inifilr: ${obj_test_inifilr}
test_orbite: ${obj_test_orbite}

depend ${makefile_dir}/depend.mk:
	makedepf90 -Wmissing -Wconfused $(addprefix -I, ${VPATH}) -nosrc $(addprefix -u , netcdf numer_rec_95 netcdf95 nr_util jumble) ${sources} >${makefile_dir}/depend.mk

clean:
	rm -f ${execut} ${objects} log

ifeq ($(findstring $(MAKECMDGOALS), clobber depend),)
include ${makefile_dir}/depend.mk
endif

-include grep.mk
-include nag_rules.mk
