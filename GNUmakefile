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
src_test_orbite = test_orbite.f orbite.f YOMCST.f unit_nml_m.f

sources := $(sort ${src_ce0l} ${src_gcm} ${src_test_ozonecm} ${src_test_inter_barxy} ${src_test_fxhyp} ${src_test_inifilr} ${src_test_orbite})

include ${general_compiler_options_dir}/settings.mk

VPATH += $(addprefix ${makefile_dir}/, $(shell cat ${makefile_dir}/directories))
cpp_macros = CPP_IIM=16,CPP_JJM=12,CPP_LLM=11
lib_list = numer_rec_95 jumble nr_util netcdf95 netcdff

# 2. Objects and executable files

obj_ce0l := $(addsuffix .o, $(basename ${src_ce0l}))
obj_gcm := $(addsuffix .o, $(basename ${src_gcm}))
obj_test_ozonecm := $(addsuffix .o, $(basename ${src_test_ozonecm}))
obj_test_inter_barxy := $(addsuffix .o, $(basename ${src_test_inter_barxy}))
obj_test_fxhyp := $(addsuffix .o, $(basename ${src_test_fxhyp}))
obj_test_inifilr := $(addsuffix .o, $(basename ${src_test_inifilr}))
obj_test_orbite := $(addsuffix .o, $(basename ${src_test_orbite}))
objects := $(addsuffix .o, $(basename ${sources}))
execut = ce0l gcm test_ozonecm test_inter_barxy test_fxhyp test_inifilr test_orbite

# 3. Rules

all: ${execut} log
ce0l: ${obj_ce0l}
gcm: ${obj_gcm}
test_ozonecm: ${obj_test_ozonecm}
test_inter_barxy: ${obj_test_inter_barxy}
test_fxhyp: ${obj_test_fxhyp}
test_inifilr: ${obj_test_inifilr}
test_orbite: ${obj_test_orbite}

depend ${makefile_dir}/depend.mk:
	makedepf90 -free -Wmissing -Wconfused $(addprefix -I, ${VPATH}) -nosrc $(addprefix -u , netcdf numer_rec_95 netcdf95 nr_util jumble) ${sources} >${makefile_dir}/depend.mk

clean:
	rm -f ${execut} ${objects} log

ifeq ($(findstring $(MAKECMDGOALS), clobber depend),)
include ${makefile_dir}/depend.mk
endif

-include grep.mk
-include nag_rules.mk
