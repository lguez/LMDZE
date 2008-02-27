# These are compiler dependent macros, meant to be included in the
# LMDZE makefile.

# For "sxf90" 2.0 on Brodie

COMPILE.f = $(FC) $(F90FLAGS) -c
FC = sxf90

# Include flags:
inc_flags = -I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I/SXlocal/pub/netCDF/3.6.1/include -I${workdir}/lib/IOIPSL_Lionel_ah -I${workdir}/lib/Numer_Rec_Lionel_f

# Other flags which do not affect run time performance:
lang_flags = -f4 -Nw -Wf "-msg b -msg d -msg o -s"

# Flags which affect run time performance:
perf_flags =

##-Cdebug -eP -eR -Pstack -Wf "-init heap=nan -init stack=nan -K a -M zdiv flovf fxovf inv setall -stmtid"
# "-M flunf" produces an error in "jacobi"

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDFLAGS =
##-Wl "-f nan"
# "-Wl,-f nan" requires "-Wf,-K a"

LDLIBS = -L${workdir}/lib/IOIPSL_Lionel_ah -lioipsl -L/SXlocal/pub/netCDF/netCDF-3.6.1/lib -lnetcdf -L${workdir}/lib/Numer_Rec_Lionel_f -lnumer_rec
