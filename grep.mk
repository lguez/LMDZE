.PHONY:  grep

grep: ${src_ce0l_only} ${common_sources}
	grep --ignore-case --files-with-matches --word-regexp dt $^
## --extended-regexp --ignore-case
