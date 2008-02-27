# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For pgf95 6.1-4

COMPILE.f = $(FC) $(F90FLAGS) -c
FC = pgf95

# Include flags:
inc_flags = -I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I${libf_dir}/filtrez -I/usr/local/netcdf-pgi/include -I/home/guez/lib/IOIPSL_Lionel/ag -I/home/guez/lib/Numer_Rec_Lionel/d

# On Duke:
## -I/usr/local/netcdf/amd64/include

# Other flags which do not affect run time performance:
lang_flags = -Mstandard -Mfree

# Flags which affect run time performance:
perf_flags = -fastsse -O3

FFLAGS = ${inc_flags} ${perf_flags}
F90FLAGS = ${inc_flags} ${lang_flags} ${perf_flags}
LDFLAGS =

LDLIBS = -L/home/guez/lib/IOIPSL_Lionel/ag -L/usr/local/netcdf-pgi/lib -L/home/guez/lib/Numer_Rec_Lionel/d -lioipsl -lnetcdf -lnumer_rec

# On Duke:
## -L/usr/local/netcdf/amd64/lib
