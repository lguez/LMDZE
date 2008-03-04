# Needs compiled modules created by NAG.

.PHONY: objects
objects: ${objects}
# (useful for NAG module builder)

nag_fcalls_options = -calledby -locate -class
nag_cross_options = -key

# NAG general options:
nag_options = -dusty ${inc_flags}

call_graph_etat0_lim: $(filter-out netcdf95.f90, ${sources_etat0_lim})
	@nag_fcalls95 ${nag_options} ${nag_fcalls_options} -listing $@ $^

call_graph_gcm: $(filter-out netcdf95.f90, ${sources_gcm})
	@nag_fcalls95 ${nag_options} ${nag_fcalls_options} -listing $@ $^

cross_ref_etat0_lim: ${sources_etat0_lim}
	@nag_xref95 ${nag_options} ${nag_cross_options} -listing $@ $^

.PHONY: clean_call_graph
clean_call_graph:
	rm -f call_graph_etat0_lim call_graph_gcm
