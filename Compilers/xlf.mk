# These are compiler dependent macros, meant to be included in the
# LMDZE makefile.

# For XL Fortran

COMPILE.f = $(FC) $(F90FLAGS) -c
FC = xlf95

# Include flags:
inc_flags = -I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I${workdir}/IOIPSL_Lionel_y ${NETCDF} -I${workdir}/Numer_Rec_Lionel

# Other flags which do not affect run time performance:
lang_flags = -qlanglvl=95pure -qnodirective -qnoescape -qsuppress=1520-050 -qwarn64

##-qattr=full -qxref=full

# Flags which affect run time performance:
perf_flags = -qcheck -qdbg -qfloat=nans -qfloat=nomaf:rndsngl:nofold -qflttrap=overflow:zerodivide:enable -qfullpath -qinitauto=7FBFFFFF -qnooptimize -qnosave -qsigtrap -qspillsize=1024

##-O3 -qnostrict -qessl -qextchk
## "-qflttrap=invalid" gives an error in "orografi.F"

FFLAGS = ${inc_flags} -qfixed ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDFLAGS = 
##-O3 -bnoquiet

LDLIBS = -L${workdir}/IOIPSL_Lionel_y -lioipsl -L${workdir}/Numer_Rec_Lionel -lnumer_rec -L/usr/local/pub/lib -lnetcdf

##-lessl
