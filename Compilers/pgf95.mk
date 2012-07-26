# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For pgf95 version 7

netcdf_inc_dir = /usr/local/netcdf-pgi/include
netcdf_lib_dir = /usr/local/netcdf-pgi/lib

numer_rec_dir = /home/guez_local/lib/Numer_Rec_Lionel/o
netcdf95_dir = /home/guez_local/lib/NetCDF95/pgf95
IOIPSL_dir = /home/guez_local/lib/IOIPSL_Lionel/ad

# Include flags:
inc_flags = $(addprefix -I, ${libf_dir} ${libf_dir}/phylmd ${netcdf_inc_dir} ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir})

# Other flags which do not affect run time performance:
lang_flags = -Mstandard -Minform=inform -Mfree -Minfo=all -Mallocatable=95 -Mbackslash

# Flags which affect run time performance:
perf_flags = -g -Kieee -Ktrap=fp -Mbounds -Mchkfpstk -Mchkptr -Mpgicoff
##-fastsse -O3
# "-Mbounds" gives an error in module "mathelp", procedure "trans_buff", 
# for a "gcm" run.

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}
LDFLAGS = -g

LDLIBS = $(addprefix -L, ${netcdf_lib_dir} ${numer_rec_dir} ${netcdf95_dir} ${IOIPSL_dir}) -lioipsl -lnetcdf95 -lnetcdf -lnumer_rec
