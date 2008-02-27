gcm_mod.txt: ${sources_gcm}
	grep --extended-regexp --ignore-case --no-filename "^ *module" $^ >$@

etat0_lim_mod.txt: ${sources_etat0_lim}
	grep --extended-regexp --ignore-case --no-filename "^ *module" $^ >$@

grep_out.txt: ${sources_etat0_lim}
	@grep --extended-regexp --files-with-matches \# $^ >$@
## --ignore-case

.PHONY: clean_grep
clean_grep: clean
	rm -f gcm_mod.txt etat0_lim_mod.txt grep_out.txt
