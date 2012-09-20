# This is a makefile for GNU `make'.

# This makefile builds LMDZE.

# Suffixes are:
# "f90" for free format, no preprocessing
# "f" for fixed format, no preprocessing

# 1. Source files

src_root = .

VPATH = ${src_root}:${src_root}/dyn3d:${src_root}/dyn3d/Vlsplt:${src_root}/filtrez:${src_root}/phylmd:${src_root}/bibio:${src_root}/phylmd/Mobidic:${src_root}/phylmd/Orography:${src_root}/phylmd/Radlwsw:${src_root}/IOIPSL:${src_root}/IOIPSL/Stringop:${src_root}/dyn3d/Read_reanalyse:${src_root}/phylmd/Thermcell:${src_root}/phylmd/CV3_routines:${src_root}/phylmd/Conflx:${src_root}/phylmd/CV_routines:${src_root}/phylmd/Interface_surf:${src_root}/dyn3d/Dissipation:${src_root}/IOIPSL/Histcom

common_sources = abort_gcm.f90 acc.f bernoui.f calendar.f90 clesphys.f90 coefils.f90 coefpoly.f comconst.f90 comdissnew.f90 comgeom.f90 comvert.f90 conf_gcm.f90 convflu.f90 convmas.f covcont.f dimens_m.f90 dimphy.f90 dimsoil.f90 diverg.f divergf.f90 diverg_gam.f divgrad2.f90 dynredem0.f90 dynredem1.f90 eigen_sort.f enercin.f ener.f90 errioipsl.f90 exner_hyb.f90 filtreg.f90 find_str.f90 flumass.f fxhyp.f fxy.f90 fxyhyper.f90 fxysinus.f fyhyp.f geopot.f90 grad.f90 gradiv2.f90 grid_change.f90 heavyside.f histbeg_totreg.f90 histclo.f90 histcom_var.f90 histdef.f90 histend.f90 histhori_regular.f90 histsync.f90 histvert.f90 indicesol.f90 iniadvtrac.f90 iniconst.f90 inidissip.f90 inifgn.f inifilr.f90 inigeom.f90 ioipslmpp.f90 jacobi.f laplacien.f90 laplacien_gam.f laplacien_rot.f laplacien_rotgam.f massbar.f massbarxy.f massdair.f mathelp.f90 minmax.f90 nxgrad.f nxgrad_gam.f nxgraro2.f90 paramet_m.f90 phyredem.f90 pressure_var.f90 q_sat.f90 rotat.f rotatf.f rotat_nfil.f scopy.f90 serre.f90 ssum.f90 strlowercase.f90 temps.f90 tourpot.f unit_nml_m.f90 vitvert.f90 ioget_calendar.f90

sources_etat0_lim = ${common_sources} caldyn0.f90 conf_dat2d.f90 conf_dat3d.f90 etat0.f90 etat0_lim.f90 extrapol.f flincom.f90 grid_atob.f90 grid_noro_m.f90 gr_int_dyn_m.f90 inter_barxy.f90 limit.f90 mva9.f90 regr_lat_time_coefoz.f90 regr_pr_o3.f90 sortvarc0.f90 startdyn.f90 start_init_orog_m.f90 start_init_phys_m.f90 start_inter_3d.f90

sources_gcm = ${common_sources} aaam_bud.f90 academic.f90 adaptdt.f addfi.f90 advect.f90 advn.f advtrac.f90 advx.f advxp.f advy.f advyp.f advz.f advzp.f aeropt.f ajsec.f90 albedo.f albsno_m.f90 bilan_dyn.f90 caladvtrac.f90 calbeta.f calcul_fluxs.f90 caldyn.f90 calfis.f90 calltherm.f90 chem.f90 clcdrag.f90 clesphys2.f90 clift.f clmain.f90 clouds_gno.f clqh.f90 cltrac.f cltracrn.f90 clvent.f90 cmpblank.f90 coefcdrag.f90 coefkz2.f90 coefkz.f90 coefkzmin.f90 comfisrtilp.f90 comgeomphy.f90 com_io_dyn.f90 concvl.f90 condsurf.f conema3.f90 conema3_m.f90 conf_guide.f90 conf_interface.f90 conflx.f90 conf_phys.f90 convect3.f coordij.f correctbid.f90 covnat.f90 ctherm.f90 cv3_closure.f90 cv3_compress.f90 cv3_feed.f90 cv3_mixing.f90 cv3_param.f90 cv3_prelim.f90 cv3_tracer.f90 cv3_trigger.f90 cv3_uncompress.f90 cv3_undilute1.f90 cv3_undilute2.f90 cv3_unsat.f90 cv3_yield.f90 cv_closure.f90 cv_compress.f90 cv_driver.f90 cv_feed.f90 cv_flag.f90 cvflag.f90 cvltr.f90 cv_mixing.f90 cvparam3.f90 cv_param.f90 cvparam.f90 cv_prelim.f90 cv_thermo.f90 cvthermo.f90 cv_trigger.f90 cv_uncompress.f90 cv_undilute1.f90 cv_undilute2.f90 cv_unsat.f90 cv_yield.f90 diagcld1.f90 diagcld2.f90 diagetpq.f90 diagphy.f90 dissip.f90 dqthermcell2.f90 dqthermcell.f90 drag_noro.f90 dteta1.f90 dudv1.f dudv2.f dump2d.f dvthermcell2.f90 dynetat0.f90 FCTTRE.f90 find_sig.f90 fisrtilp.f90 fluxstokenc.f90 flxadjtq.f90 flxasc.f90 flxbase.f90 flxddraf.f90 flxdlfs.f90 flxdtdq.f90 flxflux.f90 flxini.f90 flxmain.f90 fonte_neige.f90 formcoord.f gath_cpl.f90 gcm.f90 gensig.f90 getincom2.f90 getincom.f90 getparam.f90 gradsdef.f90 gr_fi_ecrit.f90 groupe.f90 groupeun.f gr_phy_write_3d.f90 gr_u_scal.f gr_v_scal.f guide.f90 gwprofil.f90 gwstress.f90 hbtm.f90 hgardfou.f90 histvar_seq.f90 histwrite.f90 histwrite_real.f90 inigrads.f90 ini_histday.f90 ini_histhf3d.f90 ini_histhf.f90 ini_histins.f90 ini_histrac.f90 initdynav.f90 init_dynzon.f90 initfluxsto.f90 inithist.f90 initial0.f initphysto.f90 initrrnpb.f integrd.f90 interface_surf.f90 interfoce_lim.f90 interfoce_slab.f90 interfsurf_hq.f90 interfsur_lim.f90 interpost.f interpre.f ismax.f ismin.f leapfrog.f90 lift_noro.f90 limx.f limy.f limz.f lnblnk.f lwb.f lwbv.f lwc.f lw.f lwtt.f lwttm.f lwu.f lwvb.f lwvd.f lwv.f lwvn.f mathop.f90 minmaxqfi.f90 nat2gcm.f90 newmicro.f90 nflxtr.f90 nocomma.f90 nuagecom.f90 nuage.f90 o3_chem.f90 oasis_m.f90 orbite.f90 orodrag.f90 orolift.f90 orosetup.f90 ozonecm.f90 pentes_ini.f phyetat0.f90 physiq.f90 phystokenc.f90 phytrac.f90 ppm3d.f prather.f pres2lev.f press_coefoz.f90 printflag.f90 PVtheta.f qcheck.f90 qminimum.f90 raddim.f90 raddimlw.f90 radepsi.f90 radiornpb.f90 radlwsw.f90 radopt.f90 read_reanalyse.f90 readsulfate.f reanalyse2nat.f90 regr_pr_coefoz.f90 regr_pr_comb_coefoz.f90 screenc.f90 screenp.f90 soil.f sortvarc.f90 stdlevvar.f90 sugwd.f90 suphec.f90 sw1s.f sw2s.f swclr.f swde.f sw.f90 swr.f swtt1.f swtt.f swu.f tau2alpha.f90 tetalevel.f thermcell.f90 tlift.f90 tourabs.f tracstoke.f90 transp.f transp_lay.f ustarhb.f90 vdif_kcay.f90 vlsplt.f90 vlspltqs.f90 vlx.f90 vlxqs.f90 vly.f90 vlyqs.f90 vlz.f90 wrgrads.f writedynav.f90 writehist.f yamada4.f90 yamada.f YOECUMF.f90 YOEGWD.f90 yoethf.f90 YOMCST.f90 zilch.f mathop2.f90

sources := $(sort ${sources_etat0_lim} ${sources_gcm})
# (sort so that each filename appears only once)

# 2. Objects and executable file

# Object files:
obj_etat0_lim := $(addsuffix .o, $(sort $(basename ${sources_etat0_lim})))
obj_gcm := $(addsuffix .o, $(sort $(basename ${sources_gcm})))
objects := $(addsuffix .o, $(basename ${sources}))

# Executable files:
execut = etat0_lim gcm

# 3. Compiler-dependent part:
-include ../Compilers/${FC}.mk

# 4. Rules

# Extend known suffixes:
COMPILE.f90 = $(FC) $(F90FLAGS) $(TARGET_ARCH) -c
%.o: %.f90
	$(COMPILE.f90) $(OUTPUT_OPTION) $<

.PHONY: all clean trace clobber depend

all: ${execut}

etat0_lim: ${obj_etat0_lim}
	@echo "$(FC) $(LDFLAGS) \$$^ $(LDLIBS) -o $@"
	@$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

gcm: ${obj_gcm}
	@echo "$(FC) $(LDFLAGS) \$$^ $(LDLIBS) -o $@"
	@$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@

depend ${src_root}/depend.mk:
	makedepf90 -Wmissing -Wconfused -I${VPATH} -nosrc $(addprefix -u, netcdf numer_rec_95 netcdf95 nr_util jumble) ${sources} >${src_root}/depend.mk

${src_root}/TAGS: ${sources}
	ctags -e --language-force=fortran -f $@ $^

clean:
	@echo "Removing object files."
	@rm -f ${objects}
	rm -f ${execut} trace

clobber: clean
	rm -f *.mod ${src_root}/depend.mk ${src_root}/TAGS

trace:
	echo -e "FC = ${FC}\nFFLAGS = ${FFLAGS}\nLDLIBS = ${LDLIBS}\nLDFLAGS = ${LDFLAGS}" >trace.mk

# Dependencies between object files and include files:
include ${src_root}/depend.mk

# Other rules, optionnally:
-include grep.mk
-include nag_rules.mk
