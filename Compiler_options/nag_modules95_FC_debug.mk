# These are compiler dependent macros, meant to be included in the
# makefile for LMDZE.

# For the NAGWare Fortran 95 module builder and other NAG tools

# Include flags:

FFLAGS = $(addprefix -I, ${HOME}/include ${HOME}/lib/Numer_Rec_95 ${HOME}/lib/NR_util ${HOME}/lib/NetCDF95 ${HOME}/lib/Jumble)

# NAG general options:
FFLAGS += -dusty -mismatch_all -free

nag_fcalls_options = -calledby -lines -locate -class ${FFLAGS}
nag_cross_options = -key ${FFLAGS}
version_flag = -v
