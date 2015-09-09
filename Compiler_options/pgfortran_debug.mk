# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /opt/netcdf3/pgf/include
netcdf_lib_dir = /opt/netcdf3/pgf/lib

numer_rec_95_dir = /data/guez/lib/Numer_Rec_95_pgf95_debug
nr_util_dir = /data/guez/lib/NR_util_pgf95_debug
netcdf95_dir = /data/guez/lib/NetCDF95_pgf95_debug
jumble_dir = /data/guez/lib/Jumble_pgf95_debug

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

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}) -ljumble -lnetcdf95 -lnetcdf -lnumer_rec_95 -lnr_util

version_flag = -V
