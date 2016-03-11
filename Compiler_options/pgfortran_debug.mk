# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /opt/netcdf43/pgf2013/include
netcdf_lib_dir = /opt/netcdf43/pgf2013/lib

numer_rec_95_dir = /data/guez/lib/Numer_Rec_95_pgfortran_debug
nr_util_dir = /data/guez/lib/NR_util_pgfortran_debug
netcdf95_dir = /data/guez/lib/NetCDF95_pgfortran_debug
jumble_dir = /data/guez/lib/Jumble_pgfortran_debug

lib_dir_list = ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

# Overall options:
FFLAGS += -Minform=inform

# Optimization options:
FFLAGS += -Mframe

# Debugging options:
FFLAGS += -g -gopt -Mchkfpstk -Mchkptr -Mpgicoff
# "-Mbounds" gives an error in module "mathelp", procedure "trans_buff", 
# for a "gcm" run.

# Language options:
FFLAGS += -Mfree -Mstandard -Mallocatable=95 -Mbackslash -Mdefaultunit -Mrecursive

# Target-specific Options:
FFLAGS += -Kieee -Ktrap=fp

LDFLAGS = -g

LDLIBS = $(addprefix -L, ${lib_dir_list}) -lnetcdf95 -lnetcdff -lnetcdf -lnumer_rec_95 -ljumble -lnr_util $(addprefix -rpath , ${lib_dir_list})

version_flag = -V
