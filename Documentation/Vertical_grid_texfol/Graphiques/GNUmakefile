gp_files = b_s disvert f_x strato dz_ds

odg_files = vertical_s counter_ex

compare_sampl_files = compare_sampl_9 compare_sampl_19 compare_sampl_30 compare_sampl_41

objects := $(addsuffix .tex, ${gp_files} ${compare_sampl_files}) $(addsuffix .eps, ${gp_files} ${odg_files} ${compare_sampl_files} test_disvert) $(addsuffix .pdf, ${gp_files} ${odg_files} ${compare_sampl_files} test_disvert)

all: ${objects}

%.eps %.tex: %.gp
	gnuplot $<

compare_sampl_%.eps compare_sampl_%.tex: test_disvert_%_LMD5.csv test_disvert_%_param.csv test_disvert_%_strato.csv compare_sampl.gp
	echo 'call "compare_sampl.gp" '"$*" | gnuplot -

test_disvert.eps test_disvert.pdf: half_level_6e4.csv full_level_6e4.csv half_level_7e4.csv full_level_7e4.csv
	test_disvert.py

%.pdf %.eps: %.odg
	unoconv --doctype=graphics $<
	unoconv --doctype=graphics --format=eps $<

%.pdf: %.eps
	epstopdf $<

clean:
	rm -f ${objects}
