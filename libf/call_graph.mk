# Needs compiled modules created by NAG.

.PHONY: objects
objects: ${objects}
# (useful for NAG module builder)

call_graph_etat0_lim: ${sources_etat0_lim}
	@nag_fcalls95 -calledby -dusty -locate ${inc_flags} -listing $@ $^

call_graph_gcm: $(filter-out netcdf95.f90, ${sources_gcm})
	@nag_fcalls95 -calledby -dusty -locate -class ${inc_flags} -listing $@ $^

.PHONY: clean_call_graph
clean_call_graph:
	rm -f call_graph_etat0_lim call_graph_gcm
