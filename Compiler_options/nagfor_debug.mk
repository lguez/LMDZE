# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

netcdf_dir = /opt/netcdf3/nag
numer_rec_95_dir = /data/guez/lib/Numer_Rec_95_nagfor
nr_util_dir = /data/guez/lib/NR_util_nagfor
netcdf95_dir = /data/guez/lib/NetCDF95_nagfor
jumble_dir = /data/guez/lib/Jumble_nagfor

# Include flags:
FFLAGS = $(addprefix -I, ${netcdf_dir}/include ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir})

FFLAGS += -free

LDLIBS = $(addprefix -L, ${netcdf_dir}/lib ${numer_rec_95_dir} ${netcdf95_dir} ${nr_util_dir} ${jumble_dir}) -ljumble -lnetcdf95 -lnetcdf -lnumer_rec_95 -lnr_util

version_flag = -V
