# This is the makefile for LMDZE. For GNU Make 3.80.

# Suffixes are:
# "f90" for free format, no preprocessing
# "f" for fixed format, no preprocessing

# Extend known suffixes:
%.o: %.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) $<
COMPILE.f90 = $(FC) $(F90FLAGS) -c

# Directories:
libf_dir = .

VPATH = ${libf_dir}:${libf_dir}/dyn3d:${libf_dir}/filtrez:${libf_dir}/phylmd:${libf_dir}/bibio:${libf_dir}/phylmd/Mobidic

# Compiler-dependent part:
-include compiler.mk

# Source files:
sources_etat0_lim = abort_gcm.f90 acc.f dimens_m.f90 comvert.f90 comconst.f90 paramet_m.f90 bernoui.f caldyn0.f90 coefpoly.f comdissnew.f90 conf_dat2d.f90 conf_dat3d.f90 conf_gcm.f90 convflu.f convmas.f covcont.f cray.f90 dimphy.f90 diverg.f diverg_gam.f divergf.f divgrad.f divgrad2.f dynredem0.f90 dynredem1.f90 eigen_sort.f enercin.f indicesol.f90 etat0_lim.f90 etat0.f90 exner_hyb.f90 extrapol.f filtreg.f flumass.f fxhyp.f fxy.f90 fxyhyper.f90 fxysinus.f fyhyp.f geopot.f90 grid_change.f90 grad.f gradiv.f gradiv2.f grid_atob.f90 heavyside.f iniconst.f90 inidissip.f inifgn.f inifilr.f inigeom.f90 inter_barxy.f90 jacobi.f laplacien.f laplacien_gam.f laplacien_rot.f laplacien_rotgam.f limit.f90 massbar.f massbarxy.f massdair.f minmax.f90 nxgrad.f nxgrad_gam.f nxgraro2.f nxgrarot.f phyredem.f pression.f90 psextbar.f q_sat.f90 ran1.f rotat.f rotat_nfil.f rotatf.f sortvarc0.f tourpot.f vitvert.f logic.f90 comgeom.f90 serre.f90 dimsoil.f90 temps.f90 clesphys.f90 iniprint.f90 ener.f90 advtrac_m.f90 startdyn.f90 regr1_step_av.f90 regr_lat_time_coefoz.f90 pressure_m.f90 regr1_lint.f90 start_init_phys_m.f90 start_init_orog_m.f90 grid_noro_m.f90 gr_int_dyn_m.f90 parafilt.f90 regr_pr.f90 regr3_lint.f90 regr_pr_o3.f90

sources_gcm = dimens_m.f90 comvert.f90 comconst.f90 paramet_m.f90 PVtheta.f dimphy.f90 aaam_bud.f abort_gcm.f90 acc.f adaptdt.f addfi.f advect.f advn.f advtrac.f90 advx.f advxp.f advy.f advyp.f advz.f advzp.f aeropt.f ajsec.f albedo.f atm2geo.f bernoui.f caladvtrac.f90 caldyn.f calfis.f90 calltherm.f indicesol.f90 clcdrag.f90 clift.f interface_surf.f90 clmain.f calbeta.f clqh.f clvent.f coefkz.f coefkz2.f clouds_gno.f cltrac.f cltracrn.f90 coefcdrag.f90 coefkzmin.f coefpoly.f comdissnew.f90 concvl.f condsurf.f conema3.f conema3_m.f90 conf_gcm.f90 getparam.f90 conf_guide.f conf_phys.f90 conflx.f convect3.f convflu.f convmas.f coordij.f covcont.f cray.f90 cv3_routines.f cv_driver.f cv_routines.f cvltr.f diagedyn.f diagphy.f dissip.f diverg.f diverg_gam.f divergf.f divgrad.f divgrad2.f dteta1.f dudv1.f dudv2.f dump2d.f dynetat0.f90 dynredem0.f90 dynredem1.f90 eigen_sort.f enercin.f exner_hyb.f90 filtreg.f fisrtilp.f flumass.f fluxstokenc.f formcoord.f fxhyp.f fxy.f90 fxyhyper.f90 fxysinus.f fyhyp.f gcm.f90 geopot.f90 grid_change.f90 gr_fi_ecrit.f90 grad.f gradiv.f gradiv2.f groupe.f groupeun.f guide.f90 hbtm.f heavyside.f hgardfou.f90 ini_undefSTD.f iniconst.f90 inidissip.f inifgn.f inifilr.f inigeom.f90 inigrads.f90 iniphysiq.f initdynav.f90 initfluxsto.f inithist.f90 initial0.f initphysto.f initrrnpb.f integrd.f interpost.f interpre.f ismax.f ismin.f jacobi.f laplacien.f laplacien_gam.f laplacien_rot.f laplacien_rotgam.f leapfrog.f90 limx.f limy.f limz.f lnblnk.f massbar.f massbarxy.f massdair.f minmax.f90 moy_undefSTD.f newmicro.f nflxtr.f nuage.f nxgrad.f nxgrad_gam.f nxgraro2.f nxgrarot.f orbite.f90 orografi.f ozonecm.f pentes_ini.f phyetat0.f90 phyredem.f physiq.f90 phystokenc.f phytrac.f90 plevel.f ppm3d.f prather.f pres2lev.f pression.f90 printflag.f psextbar.f q_sat.f90 qminimum.f radlwsw.f ran1.f read_reanalyse.f readsulfate.f rotat.f rotat_nfil.f rotatf.f screenc.f90 screenp.f90 soil.f sortvarc.f stdlevvar.f90 suphec.f90 tetalevel.f thermcell.f tlift.f90 tourabs.f tourpot.f transp.f transp_lay.f undefSTD.f ustarhb.f vdif_kcay.f vitvert.f vlsplt.f vlspltqs.f wrgrads.f yamada.f yamada4.f zilch.f logic.f90 comgeom.f90 serre.f90 dimsoil.f90 temps.f90 clesphys.f90 iniprint.f90 com_io_dyn.f90 tracstoke.f90 ener.f90 advtrac_m.f90 YOMCST.f90 comgeomphy.f90 ctherm.f90 writedynav.f bilan_dyn.f covnat.f gr_u_scal.f gr_v_scal.f academic.f90 oasis_m.f90 raddim.f90 radepsi.f90 radopt.f90 YOETHF.f90 FCTTRE.f90 gradsdef.f90 comfisrtilp.f90 nuagecom.f90 regr_pr_comb_coefoz.f90 regr_pr_coefoz.f90 ini_hist.f90 o3_chem.f90 comdissip.f90 conf_interface.f90 gath_cpl.f90 albsno_m.f90 parafilt.f90 regr_pr.f90

sources_test_inter_b = test_inter_barxy.f90 comconst.f90 conf_gcm.f90 iniconst.f90 inigeom.f90 iniprint.f90 logic.f90 clesphys.f90 serre.f90 comdissnew.f90 comvert.f90 comgeom.f90 fxysinus.f fxy.f90 fxyhyper.f90 cray.f90 fyhyp.f fxhyp.f inter_barxy.f90 coefpoly.f heavyside.f minmax.f90 dimens_m.f90 paramet_m.f90 temps.f90

sources_test_disv = test_disvert.f90 comvert.f90 dimens_m.f90 comconst.f90 new_unit.f90

sources := $(sort ${sources_etat0_lim} ${sources_gcm} ${sources_test_inter_b} ${sources_test_disv})
# (sort so that each filename appears only once)

# Object files:
obj_etat0_lim := $(addsuffix .o, $(sort $(basename ${sources_etat0_lim})))
obj_gcm := $(addsuffix .o, $(sort $(basename ${sources_gcm})))
obj_test_inter_b := $(addsuffix .o, $(basename ${sources_test_inter_b}))
obj_test_disv := $(sources_test_disv:.f90=.o)
objects := $(addsuffix .o, $(basename ${sources}))

# Executable files:
executables = etat0_lim gcm test_inter_barxy test_disvert

# Rules:
.PHONY: all
all: ${executables}

etat0_lim: ${obj_etat0_lim}
	@echo "$(FC) $(LDFLAGS) \$$^ $(LDLIBS) -o $@"
	@$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

gcm: ${obj_gcm}
	@echo "$(FC) $(LDFLAGS) \$$^ $(LDLIBS) -o $@"
	@$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_inter_barxy: ${obj_test_inter_b}
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

test_disvert: ${obj_test_disv}
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

TAGS: ${sources}
	etags -l fortran $^

.PHONY: clean clobber

clean:
	@echo "Removing object files."
	@rm -f ${objects}
	rm -f ${executables}

clobber: clean
	rm -f *.mod ${libf_dir}/depend.mk TAGS

.PHONY: depend
depend ${libf_dir}/depend.mk:
	makedepf90 -Wmissing -Wconfused -I ${VPATH} -u IOIPSL -u netcdf -u typesizes -u nrutil -u nrtype -u interpolation -u netcdf95 -nosrc ${sources} >${libf_dir}/depend.mk

# Dependencies between object files and include files:
include ${libf_dir}/depend.mk

# Other rules, optionnally:
-include grep.mk
-include nag_rules.mk
