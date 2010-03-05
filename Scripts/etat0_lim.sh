# This is a script in Bash.

# This script runs "etat0_lim".

trap 'echo -e \\a; exit 1' ERR
set -x

if [[ ! -f ECPHY.nc ]]
    then
    ln -sf ECDYN.nc ECPHY.nc
fi

# This is only useful if the program was compiled with G95:
G95_FPU_UNDERFLOW=${G95_FPU_UNDERFLOW:+No}

rm -f limit.nc start.nc startphy.nc coefoz_LMDZ.nc
# (If these are symbolic links then the Fortran program might not be
# able to replace them.)

time etat0_lim <etat0_lim_nml.txt >etat0_lim_out.txt 2>etat0_lim_err.txt
echo -e '\a' # beep
