# This is an extension to the LMDZE makefile, useful for NAG Fortran
# tools. It should be included in the LMDZE makefile.

# Compiled modules created by NAG are needed.

.PHONY: all_nag objects

all_nag: objects CG_etat0_lim CG_gcm CR_etat0_lim CR_gcm
objects: ${objects}

# Call graphs:
CG_etat0_lim: ${sources_etat0_lim}
	nag_fcalls95 ${nag_fcalls_options} -listing $@ $^

CG_gcm: ${sources_gcm}
	nag_fcalls95 ${nag_fcalls_options} -listing $@ $^

# Cross references:
CR_etat0_lim: ${sources_etat0_lim}
	nag_xref95 ${nag_cross_options} -listing $@ $^

CR_gcm: ${sources_gcm}
	nag_xref95 ${nag_cross_options} -listing $@ $^

CG_etat0_lim CG_gcm CR_etat0_lim CR_gcm: objects

.PHONY: clean_nag
clean_nag:
	rm -f CG_etat0_lim CG_gcm CR_etat0_lim CR_gcm
