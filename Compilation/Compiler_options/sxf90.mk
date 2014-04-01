# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For "sxf90" 2.0

netcdf_inc_dir = /SXlocal/pub/netCDF/3.6.1/include
netcdf_lib_dir = /SXlocal/pub/netCDF/3.6.1/lib

numer_rec_dir = ${workdir}/lib/Numer_Rec_Lionel/x
netcdf95_dir = ${workdir}/lib/NetCDF95_sxf90
IOIPSL_dir = ${workdir}/lib/IOIPSL_Lionel_am

# Include flags:
inc_flags = $(addprefix -I, ${libf_dir} ${libf_dir}/dyn3d ${libf_dir}/phylmd ${libf_dir}/filtrez ${netcdf_inc_dir} ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir})

# Other flags which do not affect run time performance:
lang_flags = -f4 -Nw -Wf "-msg b -msg d -msg o -s"

# Flags which affect run time performance:
perf_flags = -Cdebug -eP -eR -Pstack -Wf "-init heap=nan -init stack=nan -K a -M zdiv flovf fxovf inv setall -stmtid"

# "-M flunf" produces an error in "jacobi"

FFLAGS = ${inc_flags} ${lang_flags} ${perf_flags}

LDFLAGS = -Wl "-f nan"
# "-Wl,-f nan" requires "-Wf,-K a"

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir}) -lioipsl -lnetcdf95 -lnetcdf -lnumer_rec
