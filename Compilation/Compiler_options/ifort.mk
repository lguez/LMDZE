# These are machine dependent macros, meant to be included in the
# LMDZE makefile

# For ifort 9

FFLAGS = -I${libf_dir} -I${libf_dir}/dyn3d -I${libf_dir}/phylmd -I/home/guez/NetCDF/netcdf-3.6.1_ifort/include -I/home/guez/lib/IOIPSL/IOIPSL_t -I/home/guez/lib/Numerical_Recipes_Lionel/ifort -assume minus0 -check all -debug extended -debug-parameters all -error_limit 1 -fltconsistency -fpe0 -fpstkchk -ftrapuv -g -inline-debug-info -O0 -traceback -free -warn all -warn stderrors

LDFLAGS = 

LDLIBS = -L/home/guez/NetCDF/netcdf-3.6.1_ifort/lib -L/home/guez/lib/IOIPSL/IOIPSL_t -L/home/guez/lib/Numerical_Recipes_Lionel/ifort -lioipsl -lnetcdf -lnumer_rec
