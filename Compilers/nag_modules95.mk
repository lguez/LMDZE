# These are machine dependent macros, meant to be included in the
# LMDZE makefile.

# For the NAGWare Fortran 95 module builder

# Include flags:
inc_flags = -I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I/home/guez/NetCDF/netcdf-3.6.1_NAG -I/home/guez/lib/IOIPSL_Lionel/ai -I/home/guez/lib/Numer_Rec_Lionel/g

COMPILE.f90 = nag_modules95.sh $@ ${inc_flags} -dusty -mismatch_all
COMPILE.f = ${COMPILE.f90}
OUTPUT_OPTION =
