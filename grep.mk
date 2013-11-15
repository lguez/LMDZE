gcm_mod.txt: ${sources_gcm}
	grep --extended-regexp --ignore-case --no-filename "^ *module" $^ >$@

ce0l_mod.txt: ${sources_ce0l}
	grep --extended-regexp --ignore-case --no-filename "^ *module" $^ >$@

grep: ${sources_ce0l}
	@echo grep
	@grep --ignore-case --files-with-matches --word-regexp dt $^
## --extended-regexp --ignore-case

.PHONY: clean_grep grep
clean_grep: clean
	rm -f gcm_mod.txt ce0l_mod.txt
