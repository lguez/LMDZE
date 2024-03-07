sources = leapfrog time_gcm masque sollw ECDYN phis Traceurs tourpot ps_gcm Grilles_horizontales_LMDZ bucket ECPHY qsol day_ini xprimu jfilt rlon_rlat dyn_phy icb iflag sigt plcl pctsrf temperature tsol inter_barxy lecture_fichiers_ce0l vorticity circulation thermcell wa surface_energy zvert

objects := $(addsuffix .pdf, ${sources})

all: ${objects}

%.pdf: %.odg
	unoconv --doctype=graphics $<

%.pdf: %.gv
	dot -Tpdf -o $@ $<

%.pdf: %.py
	$<

clean:
	rm -f ${objects}
