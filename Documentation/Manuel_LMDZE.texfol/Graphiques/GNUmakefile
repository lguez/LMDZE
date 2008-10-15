gp_files = b_s disvert f_x strato1 strato2

odg_files = Cellule leapfrog libraries ozone regridding_proced regr_latit regr_lat_time_coefoz regr_pr_coefoz regr_p_stepav regr_p_lint Sigma_interp_vert

compare_sampl_files = compare_sampl_9 compare_sampl_19 compare_sampl_30 compare_sampl_41

objects := $(addsuffix .tex, ${gp_files} ${compare_sampl_files}) $(addsuffix .eps, ${gp_files} ${odg_files} ${compare_sampl_files}) $(addsuffix .pdf, ${gp_files} ${odg_files} ${compare_sampl_files})

%.eps %.tex: %.gp
	gnuplot $<
	sed -i 's/\includegraphics{/\includegraphics{Graphiques\//' $*.tex

compare_sampl_%.eps compare_sampl_%.tex: test_disvert_%_LMD5.csv test_disvert_%_param.csv test_disvert_%_strato1.csv test_disvert_%_strato2.csv compare_sampl.gp
	echo 'call "compare_sampl.gp" '"$*" | gnuplot -
	sed -i 's/\includegraphics{/\includegraphics{Graphiques\//' compare_sampl_$*.tex

# Join rules for pdf and eps from odg to avoid circular dependency
# between eps and pdf:
%.pdf %.eps: %.odg
	unoconv $<
	convert $*.pdf $*.eps

%.pdf: %.eps
	epstopdf $<

all: ${objects}

clean:
	rm -f ${objects}
