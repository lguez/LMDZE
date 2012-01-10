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

VPATH = ${libf_dir}:${libf_dir}/dyn3d:${libf_dir}/dyn3d/Vlsplt:${libf_dir}/filtrez:${libf_dir}/phylmd:${libf_dir}/bibio:${libf_dir}/phylmd/Mobidic:${libf_dir}/phylmd/Orography:${libf_dir}/phylmd/Radlwsw:${libf_dir}/IOIPSL:${libf_dir}/IOIPSL/Stringop:${libf_dir}/dyn3d/Read_reanalyse:${libf_dir}/phylmd/Thermcell:${libf_dir}/phylmd/CV3_routines:${libf_dir}/phylmd/Conflx:${libf_dir}/phylmd/CV_routines:${libf_dir}/phylmd/Interface_surf:${libf_dir}/dyn3d/Dissipation

# Compiler-dependent part:
-include compiler.mk

# Source files:
sources_etat0_lim = acc.f dimens_m.f90 comvert.f90 comconst.f90 paramet_m.f90 bernoui.f caldyn0.f90 coefpoly.f comdissnew.f90 conf_dat2d.f90 conf_dat3d.f90 conf_gcm.f90 convflu.f convmas.f covcont.f dimphy.f90 diverg.f diverg_gam.f divergf.f divgrad.f divgrad2.f dynredem0.f90 dynredem1.f90 eigen_sort.f enercin.f indicesol.f90 etat0_lim.f90 etat0.f90 exner_hyb.f90 extrapol.f filtreg.f90 flumass.f fxhyp.f fxy.f90 fxyhyper.f90 fxysinus.f fyhyp.f geopot.f90 grid_change.f90 grad.f gradiv.f gradiv2.f90 grid_atob.f90 heavyside.f iniconst.f90 inidissip.f90 inifgn.f inifilr.f90 inigeom.f90 inter_barxy.f90 jacobi.f laplacien.f laplacien_gam.f laplacien_rot.f laplacien_rotgam.f limit.f90 massbar.f massbarxy.f massdair.f minmax.f90 nxgrad.f nxgrad_gam.f nxgraro2.f nxgrarot.f phyredem.f90 q_sat.f90 rotat.f rotat_nfil.f rotatf.f sortvarc0.f90 tourpot.f vitvert.f90 logic.f90 comgeom.f90 serre.f90 dimsoil.f90 temps.f90 clesphys.f90 iniprint.f90 ener.f90 iniadvtrac.f90 startdyn.f90 regr_lat_time_coefoz.f90 pressure_var.f90 start_init_phys_m.f90 start_init_orog_m.f90 grid_noro_m.f90 gr_int_dyn_m.f90 regr_pr_o3.f90 coefils.f90 scopy.f90 ssum.f90 histcom.f90 ioipslmpp.f90 errioipsl.f90 histcom_var.f90 mathelp.f90 calendar.f90 flincom.f90 find_str.f90 strlowercase.f90 abort_gcm.f90 mva9.f90 start_inter_3d.f90

sources_gcm = dimens_m.f90 comvert.f90 comconst.f90 paramet_m.f90 PVtheta.f dimphy.f90 aaam_bud.f90 abort_gcm.f90 acc.f adaptdt.f addfi.f90 advect.f90 advn.f advtrac.f90 advx.f advxp.f advy.f advyp.f advz.f advzp.f aeropt.f ajsec.f90 albedo.f bernoui.f caladvtrac.f90 caldyn.f90 calfis.f90 calltherm.f90 indicesol.f90 clcdrag.f90 clift.f interface_surf.f90 clmain.f90 calbeta.f clqh.f90 clvent.f coefkz.f90 coefkz2.f90 clouds_gno.f cltrac.f cltracrn.f90 coefcdrag.f90 coefkzmin.f90 coefpoly.f comdissnew.f90 concvl.f90 condsurf.f conema3.f90 conema3_m.f90 conf_gcm.f90 getparam.f90 conf_guide.f90 conf_phys.f90 conflx.f90 convect3.f convflu.f convmas.f coordij.f covcont.f cv_driver.f90 cvltr.f diagphy.f dissip.f90 diverg.f diverg_gam.f divergf.f divgrad.f divgrad2.f dteta1.f90 dudv1.f dudv2.f dump2d.f dynetat0.f90 dynredem0.f90 dynredem1.f90 eigen_sort.f enercin.f exner_hyb.f90 filtreg.f90 fisrtilp.f90 flumass.f fluxstokenc.f90 formcoord.f fxhyp.f fxy.f90 fxyhyper.f90 fxysinus.f fyhyp.f gcm.f90 geopot.f90 grid_change.f90 gr_fi_ecrit.f90 grad.f gradiv.f gradiv2.f90 groupe.f90 groupeun.f guide.f90 hbtm.f90 heavyside.f hgardfou.f90 iniconst.f90 inidissip.f90 inifgn.f inifilr.f90 inigeom.f90 inigrads.f90 initdynav.f90 initfluxsto.f90 inithist.f90 initial0.f initphysto.f90 initrrnpb.f integrd.f90 interpost.f interpre.f ismax.f ismin.f jacobi.f laplacien.f laplacien_gam.f laplacien_rot.f laplacien_rotgam.f leapfrog.f90 limx.f limy.f limz.f lnblnk.f massbar.f massbarxy.f massdair.f minmax.f90 newmicro.f90 nflxtr.f nuage.f90 nxgrad.f nxgrad_gam.f nxgraro2.f nxgrarot.f orbite.f90 ozonecm.f90 pentes_ini.f phyetat0.f90 phyredem.f90 physiq.f90 phystokenc.f90 phytrac.f90 ppm3d.f prather.f pres2lev.f printflag.f90 q_sat.f90 qminimum.f90 read_reanalyse.f90 readsulfate.f rotat.f rotat_nfil.f rotatf.f screenc.f90 screenp.f90 soil.f sortvarc.f90 stdlevvar.f90 suphec.f90 tetalevel.f thermcell.f90 tlift.f90 tourabs.f tourpot.f transp.f transp_lay.f ustarhb.f vdif_kcay.f vitvert.f90 vlsplt.f90 vlx.f90 vly.f90 vlz.f90 vlspltqs.f90 wrgrads.f yamada.f yamada4.f90 zilch.f logic.f90 comgeom.f90 serre.f90 dimsoil.f90 temps.f90 clesphys.f90 iniprint.f90 com_io_dyn.f90 tracstoke.f90 ener.f90 iniadvtrac.f90 YOMCST.f90 comgeomphy.f90 ctherm.f90 writedynav.f90 bilan_dyn.f90 covnat.f90 gr_u_scal.f gr_v_scal.f academic.f90 oasis_m.f90 raddim.f90 radepsi.f90 radopt.f90 yoethf.f90 FCTTRE.f90 gradsdef.f90 comfisrtilp.f90 nuagecom.f90 regr_pr_comb_coefoz.f90 regr_pr_coefoz.f90 ini_histday.f90 ini_histhf3d.f90 ini_histrac.f90 ini_histhf.f90  ini_histins.f90 o3_chem.f90 conf_interface.f90 gath_cpl.f90 albsno_m.f90 YOECUMF.f90 YOEGWD.f90 raddimlw.f90 pressure_var.f90 clesphys2.f90 qcheck.f90 radiornpb.f90 minmaxqfi.f90 press_coefoz.f90 drag_noro.f90 gwstress.f90 orodrag.f90 orosetup.f90 gwprofil.f90 lift_noro.f90 orolift.f90 sugwd.f90 lw.f sw.f90 lwc.f lwu.f lwvd.f swclr.f swtt.f lwb.f lwtt.f lwv.f lwvn.f sw1s.f swde.f swtt1.f lwbv.f lwttm.f lwvb.f radlwsw.f90 sw2s.f swr.f swu.f coefils.f90 scopy.f90 ssum.f90 histcom.f90 getincom.f90 histwrite.f90 calendar.f90 errioipsl.f90 histcom_var.f90 ioipslmpp.f90 mathelp.f90 cmpblank.f90 find_str.f90 nocomma.f90 find_sig.f90 gensig.f90 strlowercase.f90 gr_phy_write_3d.f90 mathop.f90 tau2alpha.f90 correctbid.f90 nat2gcm.f90 reanalyse2nat.f90 vlxqs.f90 vlyqs.f90 histwrite_real.f90 histvar_seq.f90 diagetpq.f90 dqthermcell2.f90 dqthermcell.f90 dvthermcell2.f90 cv3_closure.f90 cv3_param.f90 cv3_uncompress.f90 cv3_yield.f90 cv3_compress.f90 cv3_prelim.f90 cv3_undilute1.f90 cv3_feed.f90 cv3_tracer.f90 cv3_undilute2.f90 cv3_mixing.f90 cv3_trigger.f90 cv3_unsat.f90 cvparam3.f90 cvthermo.f90 cvflag.f90 cv_thermo.f90 cv_flag.f90 getincom2.f90 diagcld1.f90 diagcld2.f90 chem.f90 flxasc.f90 flxddraf.f90 flxdtdq.f90 flxini.f90 flxsetup.f90 flxadjtq.f90 flxbase.f90 flxdlfs.f90 flxflux.f90 flxmain.f90 cv_closure.f90 cv_mixing.f90 cv_trigger.f90 cv_undilute2.f90 cv_compress.f90 cv_param.f90 cv_uncompress.f90 cv_unsat.f90 cv_feed.f90 cv_prelim.f90 cv_undilute1.f90 cv_yield.f90 cvparam.f90 calcul_fluxs.f90 interfoce_slab.f90 interfsur_lim.f90 fonte_neige.f90 interfoce_lim.f90 interfsurf_hq.f90 writehist.f

sources := $(sort ${sources_etat0_lim} ${sources_gcm})
# (sort so that each filename appears only once)

# Object files:
obj_etat0_lim := $(addsuffix .o, $(sort $(basename ${sources_etat0_lim})))
obj_gcm := $(addsuffix .o, $(sort $(basename ${sources_gcm})))
objects := $(addsuffix .o, $(basename ${sources}))

# Executable files:
executables = etat0_lim gcm

# Rules:
.PHONY: all
all: ${executables}

etat0_lim: ${obj_etat0_lim}
	@echo "$(FC) $(LDFLAGS) \$$^ $(LDLIBS) -o $@"
	@$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

gcm: ${obj_gcm}
	@echo "$(FC) $(LDFLAGS) \$$^ $(LDLIBS) -o $@"
	@$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

TAGS: ${sources}
	etags --language-force=fortran $^

.PHONY: clean clobber

clean:
	@echo "Removing object files."
	@rm -f ${objects}
	rm -f ${executables}

clobber: clean
	rm -f *.mod ${libf_dir}/depend.mk TAGS

.PHONY: depend
depend ${libf_dir}/depend.mk:
	makedepf90 -Wmissing -Wconfused -I ${VPATH} -u netcdf -u numer_rec -u netcdf95 -u nr_util -u jumble -nosrc ${sources} >${libf_dir}/depend.mk

# Dependencies between object files and include files:
include ${libf_dir}/depend.mk

# Other rules, optionnally:
-include grep.mk
-include nag_rules.mk
