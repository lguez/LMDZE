# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /opt/netcdf3/pgf/include
netcdf_lib_dir = /opt/netcdf3/pgf/lib

numer_rec_95_dir = /data/guez/lib/Numer_Rec_95_pgfortran
nr_util_dir = /data/guez/lib/NR_util_pgfortran
netcdf95_dir = /data/guez/lib/NetCDF95_pgfortran
jumble_dir = /data/guez/lib/Jumble_pgfortran

lib_dir_list = ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Optimization options:
FFLAGS += -Mipa=cg,required,safeall

# Language options:
FFLAGS += -Mfree -Mstandard -Mallocatable=95 -Mbackslash -Mdefaultunit -Mrecursive

LDLIBS = $(addprefix -L, ${lib_dir_list}) -lnetcdf95 -lnetcdf -lnumer_rec_95 -ljumble -lnr_util $(addprefix -rpath , ${lib_dir_list})

version_flag = -V

LDFLAGS = -Mipa=cg,required,safeall
