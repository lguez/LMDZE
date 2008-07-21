# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For the NAGWare Fortran 95 module builder and other NAG tools.

# Include flags:
inc_flags = -I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I/home/guez_local/include -I/home/guez_local/include/NetCDF_nag_modules95 -I/home/guez_local/lib/IOIPSL_Lionel/ai -I/home/guez_local/lib/Numer_Rec_Lionel/g -I/home/guez_local/lib/NetCDF95/nag_modules95

# NAG general options:
nag_gl_options = ${inc_flags} -dusty -mismatch_all

COMPILE.f90 = nag_modules95.sh $@ ${nag_gl_options}
COMPILE.f = ${COMPILE.f90}
OUTPUT_OPTION =

nag_fcalls_options = -calledby -lines -locate -class ${nag_gl_options}
nag_cross_options = -key ${nag_gl_options}
