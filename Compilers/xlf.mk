# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For IBM XL Fortran

FC = xlf95

numer_rec_dir =  ${workdir}/lib/Numer_Rec_Lionel_v
netcdf95_dir = ${workdir}/lib/NetCDF95
IOIPSL_dir = ${workdir}/lib/IOIPSL_Lionel_ar

# Include flags:
inc_flags = $(addprefix -I, ${libf_dir} ${libf_dir}/dyn3d ${libf_dir}/phylmd ${libf_dir}/filtrez ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir}) ${NETCDF_FFLAGS}

# Other flags which do not affect run time performance:
lang_flags = -qlanglvl=95pure -qnodirective -qnoescape -qsuppress=1520-050 -qwarn64

##-qattr=full -qxref=full

# Flags which affect run time performance:
perf_flags = -qdbg -qfloat=nans -qfloat=nomaf:rndsngl:nofold -qflttrap=overflow:zerodivide:enable -qfullpath -qinitauto=7FBFFFFF -qnooptimize -qnosave -qsigtrap -qspillsize=1024

##-O3 -qnostrict -qessl -qextchk
## "-qflttrap=invalid" gives an error in "orografi.F"
##-qcheck severe error in etat0
# "-qcheck -qextchk" give an error in module "mathelp", procedure
# "trans_buff", for a "gcm" run.

FFLAGS = ${inc_flags} -qfixed ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDFLAGS = 
##-O3 -bnoquiet

LDLIBS = $(addprefix -L, ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir}) -lioipsl -lnetcdf95 -lnumer_rec ${NETCDF_LDFLAGS}

##-lessl
