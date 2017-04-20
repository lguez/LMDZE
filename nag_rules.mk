# This is an extension to the LMDZE makefile, useful for NAG Fortran
# tools. It should be included in the LMDZE makefile.

# Compiled modules created by NAG are needed.

.PHONY: all_nag
all_nag: CG_ce0l CG_gcm

# Call graphs:
CG_ce0l: ${src_ce0l}
	nagfor =callgraph ${nag_fcalls_options} -o $@ $^

CG_gcm: ${src_gcm}
	nagfor =callgraph ${nag_fcalls_options} -o $@ $^

.PHONY: clean_nag
clean_nag:
	rm -f CG_ce0l CG_gcm
