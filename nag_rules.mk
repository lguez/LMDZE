# This is an extension to the LMDZE makefile, useful for NAG Fortran
# tools. It should be included in the LMDZE makefile.

# Compiled modules created by NAG are needed.

.PHONY: all_nag objects

all_nag: objects CG_ce0l CG_gcm CR_ce0l CR_gcm
objects: ${objects}

# Call graphs:
CG_ce0l: ${src_ce0l}
	nagfor =callgraph ${nag_fcalls_options} -o $@ $^

CG_gcm: ${src_gcm}
	nagfor =callgraph ${nag_fcalls_options} -o $@ $^

# Cross references:
CR_ce0l: ${src_ce0l}
	nag_xref95 ${nag_cross_options} -listing $@ $^

CR_gcm: ${src_gcm}
	nag_xref95 ${nag_cross_options} -listing $@ $^

.PHONY: clean_nag
clean_nag:
	rm -f CG_ce0l CG_gcm CR_ce0l CR_gcm
