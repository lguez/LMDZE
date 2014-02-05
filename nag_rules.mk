# This is an extension to the LMDZE makefile, useful for NAG Fortran
# tools. It should be included in the LMDZE makefile.

# Compiled modules created by NAG are needed.

.PHONY: all_nag objects

all_nag: objects CG_ce0l CG_gcm CR_ce0l CR_gcm
objects: ${objects}

# Call graphs:
CG_ce0l: ${common_sources} ${src_no_main_ce0l_only} ce0l.f90
	nag_fcalls95 ${nag_fcalls_options} -listing $@ $^

CG_gcm: ${common_sources} ${src_no_main_gcm_only} gcm.f90
	nag_fcalls95 ${nag_fcalls_options} -listing $@ $^

# Cross references:
CR_ce0l: ${common_sources} ${src_no_main_ce0l_only} ce0l.f90
	nag_xref95 ${nag_cross_options} -listing $@ $^

CR_gcm: ${common_sources} ${src_no_main_gcm_only} gcm.f90
	nag_xref95 ${nag_cross_options} -listing $@ $^

.PHONY: clean_nag
clean_nag:
	rm -f CG_ce0l CG_gcm CR_ce0l CR_gcm
