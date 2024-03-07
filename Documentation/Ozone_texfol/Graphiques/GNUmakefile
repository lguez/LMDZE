odg_files = Cellule ozone_flow regridding_proced regr_latit regr_lat_time_coefoz regr_pr_coefoz regr_p_stepav regr_p_lint Sigma_interp_vert

objects := $(addsuffix .pdf, ${odg_files})

all: ${objects}

%.pdf: %.odg
	unoconv --doctype=graphics $<

clean:
	rm -f ${objects}
