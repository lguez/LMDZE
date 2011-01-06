# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For the NAGWare Fortran 95 module builder and other NAG tools

# Include flags:
inc_flags = $(addprefix -I, ${libf_dir} ${libf_dir}/dyn3d ${libf_dir}/phylmd ${libf_dir}/filtrez /net/vierne/user/guez_local/include ${HOME}/include ${HOME}/lib/Numer_Rec_Lionel_s ${HOME}/lib/NR_util_g ${HOME}/lib/NetCDF95)

# NAG general options:
nag_gl_options = ${inc_flags} -dusty -mismatch_all

COMPILE.f90 = nag_modules95.sh $@ ${nag_gl_options}
COMPILE.f = ${COMPILE.f90}
OUTPUT_OPTION =

nag_fcalls_options = -calledby -lines -locate -class ${nag_gl_options}
nag_cross_options = -key ${nag_gl_options}
