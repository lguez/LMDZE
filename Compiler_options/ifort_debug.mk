# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For ifort 12

netcdf_dir = /smplocal/pub/NetCDF/4.1.3/seq
##netcdf_dir = /opt/netcdf42/ifort

netcdf_inc_dir = ${netcdf_dir}/include
netcdf_lib_dir = ${netcdf_dir}/lib

numer_rec_95_dir = ${workdir}/lib/Numer_Rec_95_debug
##numer_rec_95_dir = /data/guez/lib/Numer_Rec_95_ifort_debug

nr_util_dir = ${workdir}/lib/NR_util_debug
##nr_util_dir = /data/guez/lib/NR_util_ifort_debug

netcdf95_dir = ${workdir}/lib/NetCDF95_debug
##netcdf95_dir = /data/guez/lib/NetCDF95_ifort_debug

jumble_dir =${workdir}/lib/Jumble_debug
##jumble_dir =/data/guez/lib/Jumble_ifort_debug

lib_dir_list = ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Optimization:
FFLAGS += -O0

# Floating Point:
FFLAGS += -fp-stack-check -fpe-all=0

# Debug:
FFLAGS += -debug -debug-parameters all -ftrapuv

# Language:
FFLAGS += -free -noaltparam -assume minus0,noold_xor -check bounds,format,output_conversion,pointers,uninit -stand f95

# Data:
FFLAGS += -auto

# Compiler Diagnostics:
FFLAGS += -warn declarations,interfaces,stderrors,truncated_source,uncalled,unused -traceback -diag-error-limit 1

comma = ,

LDLIBS = $(addprefix -L, ${lib_dir_list}) -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec_95 -ljumble -lnr_util $(addprefix -Wl${comma}-rpath${comma}, ${lib_dir_list})

version_flag = -logo
