# This is a script in Bash.

# This script collects input files necessary for the program "gcm" and
# runs "gcm". Run this script in the directory where you want "gcm" to run.
# The "ln" command does not report missing targets, so we check them first.

trap 'echo -e \\a; exit 1' ERR

# Specify directories:

REL_dir=/home/guez/In_transit/LMDZE_work/Results_etat0_lim
##REL_dir=$workdir/LMDZE/Results_etat0_lim

IGCM_dir=~/Documents/Utilisation_LMDZ/Input_gcm
##IGCM_dir=~

executable_dir=/home/guez/In_transit/LMDZE_work/Compil_prod
##executable_dir=$workdir/LMDZE/Compil_prod

RGCM_dir=/home/guez/In_transit/LMDZE_work/Results_gcm
##RGCM_dir=$workdir/LMDZE/Results_gcm
# (used only for a restart or a comparison)

set -x

test -d $REL_dir
test -d $IGCM_dir
test -d $executable_dir

set +x

read -p \
    "Identifier for the set of parameters (\"*.def\" files and namelists)? " \
    igcm_id
read -p \
    "Run number of \"etat0_lim\", for \"limit.nc\" and \"coefoz_LMDZ.nc\"? " \
    numb_cr

read -p \
"Do you want to restart from the end of a previous \"gcm\" run? (y/[n]) " \
    restart
if [[ $restart = y ]]
    then
    echo "Previous run number of \"gcm\" for \"start.nc\" and \"startphy.nc\ ?"
    read prev_gcm
fi

set -x

my_host=`hostname`

if [[ $my_host = vierne ]]
    then
    # This is only useful if the program was compiled with G95:
    G95_FPU_UNDERFLOW=${G95_FPU_UNDERFLOW:+No}
    G95_FPU_INVALID=${G95_FPU_INVALID:+No}
fi

test -f $REL_dir/$numb_cr/limit.nc
##test -f $REL_dir/$numb_cr/coefoz_LMDZ.nc
ln -f -s $REL_dir/$numb_cr/limit.nc $REL_dir/$numb_cr/coefoz_LMDZ.nc .

if [[ $restart != y ]]
    then
    # Start-up files come from "etat0_lim":
    ln -f -s $REL_dir/$numb_cr/start*.nc .
else
    ln -f -s $RGCM_dir/$prev_gcm/restart.nc start.nc
    ln -f -s $RGCM_dir/$prev_gcm/restartphy.nc startphy.nc
fi

rm -f *.def
test -f $IGCM_dir/$igcm_id/run.def
ln -s $IGCM_dir/$igcm_id/*.def .
date

if [[ $my_host = brodie ]]
    then
    rsh brodie01 \
    export F_PROGINF=YES\; \
    cd $PWD\; \
    $executable_dir/gcm <$IGCM_dir/$igcm_id/gcm_nml.txt >gcm_out.txt \
    2>gcm_err.txt
elif [[ $my_host = zahir* ]]
    then
    hpmcount $executable_dir/gcm \
	<$IGCM_dir/$igcm_id/gcm_nml.txt >gcm_out.txt 2>gcm_err.txt
else
    time $executable_dir/gcm \
	<$IGCM_dir/$igcm_id/gcm_nml.txt >gcm_out.txt 2>gcm_err.txt
fi

echo -e '\a' # beep
set +x
read -p "Previous run number for comparison [none] ? "
if [[ -n $REPLY ]]
    then
    selective_diff.sh $RGCM_dir/$REPLY .
fi
