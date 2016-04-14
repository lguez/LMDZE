.PHONY: grep

grep: ${src_ce0l}
	grep --ignore-case --files-with-matches --word-regexp dt $^
## --extended-regexp --ignore-case

src_with_dir: ${sources}
	echo $^ >$@
