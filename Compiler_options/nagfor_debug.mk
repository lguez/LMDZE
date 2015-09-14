# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_inc_dir = /opt/netcdf3/nag/include
netcdf_lib_dir = /opt/netcdf3/nag/lib

numer_rec_95_dir = /data/guez/lib/Numer_Rec_95_nagfor
nr_util_dir = /data/guez/lib/NR_util_nagfor
netcdf95_dir = /data/guez/lib/NetCDF95_nagfor
jumble_dir = /data/guez/lib/Jumble_nagfor

lib_dir_list = ${netcdf_lib_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_inc_dir} ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

FFLAGS += -free

LDLIBS = $(addprefix -L, ${lib_dir_list}) -lnetcdf95 -lnetcdf -lnumer_rec_95 -ljumble -lnr_util

version_flag = -V
nag_fcalls_options = -calledby -index ${FFLAGS}
nag_cross_options = -key ${FFLAGS}
