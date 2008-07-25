gp_files = b_s disvert f_x strato1 strato2

odg_files = Cellule leapfrog libraries ozone regridding_proced regr_latit regr_lat_time_coefoz regr_pr_coefoz regr_p_stepav regr_p_lint

objects := $(addsuffix .tex, ${gp_files}) $(addsuffix .eps, ${gp_files} ${odg_files}) $(addsuffix .pdf, ${gp_files} ${odg_files}) compare_sampl_9.pdf compare_sampl_19.pdf compare_sampl_30.pdf compare_sampl_41.pdf

%.eps %.tex: %.gp
	gnuplot $<
	sed -i 's/\includegraphics{/\includegraphics{Graphiques\//' $*.tex

# Join rules for pdf and eps from odg to avoid circular dependancy
# between eps and pdf:
%.pdf %.eps: %.odg
	unoconv $<
	convert $*.pdf $*.eps

%.pdf: %.eps
	epstopdf $<

all: ${objects}

clean:
	rm -f ${objects}
