# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For the NAGWare Fortran 95 module builder and other NAG tools

# Include flags:
inc_flags = $(addprefix -I, ${HOME}/include ${HOME}/lib/Numer_Rec_95 ${HOME}/lib/NR_util_g ${HOME}/lib/NetCDF95 ${HOME}/lib/Jumble)

# NAG general options:
nag_gl_options = ${inc_flags} -dusty -mismatch_all

FFLAGS= ${nag_gl_options}
F90FLAGS = ${FFLAGS}

nag_fcalls_options = -calledby -lines -locate -class ${nag_gl_options}
nag_cross_options = -key ${nag_gl_options}
