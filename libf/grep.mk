gcm_mod.txt: ${sources_gcm}
	grep --extended-regexp --ignore-case --no-filename "^ *module" $^ >$@

etat0_lim_mod.txt: ${sources_etat0_lim}
	grep --extended-regexp --ignore-case --no-filename "^ *module" $^ >$@

grep: ${sources_etat0_lim}
	@echo grep
	@grep --ignore-case --files-with-matches dtvr $^
## --extended-regexp --ignore-case

.PHONY: clean_grep grep
clean_grep: clean
	rm -f gcm_mod.txt etat0_lim_mod.txt
