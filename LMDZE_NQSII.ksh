# This is a script in KornShell with directives for NQSII.

# This script collects input files necessary for the program "gcm" and
# runs "gcm". Submit this script to NQSII on Brodie.

#PBS -S /bin/ksh
#PBS -N gcm_97
#PBS -o gcm_97_out
#PBS -j o
#PBS -l cputim_job=30:00
#PBS -l memsz_job=300MB

set -x
trap 'exit 1' ERR

cd ${TMPDIR:?}

# Specify directories:
REL_dir=$workdir/LMDZE/Results_etat0_lim
IGCM_dir=~
executable_dir=$workdir/LMDZE/Compil_prod

# Identifier for the set of parameters ("*.def" files and namelists):
igcm_id=igcm34

# Run number of "etat0_lim", for "limit.nc" and "coefoz_LMDZ.nc":
numb_cr=53

# Restart from the end of a previous "gcm" run?
restart=n

if [[ $restart = y ]]
    then
    RGCM_dir=~/Documents/Utilisation_LMDZ/Results_gcm

    # Previous run number of "gcm" for "start.nc" and "startphy.nc":
    prev_gcm=
fi

test -f $REL_dir/$numb_cr/limit.nc
test -f $REL_dir/$numb_cr/coefoz_LMDZ.nc
cp -p $REL_dir/$numb_cr/limit.nc $REL_dir/$numb_cr/coefoz_LMDZ.nc .

if [[ $restart != y ]]
    then
    # Start-up files come from "etat0_lim":
    cp -p $REL_dir/$numb_cr/start*.nc .
else
    cp -p $RGCM_dir/$prev_gcm/restart.nc start.nc
    cp -p $RGCM_dir/$prev_gcm/restartphy.nc startphy.nc
fi

test -f $IGCM_dir/$igcm_id/run.def
cp -p $IGCM_dir/$igcm_id/*.def .

cp -p $executable_dir/gcm $IGCM_dir/$igcm_id/gcm_nml.txt .
ls -l
export F_PROGINF=YES
trap - ERR
gcm <gcm_nml.txt
## >gcm_out.txt
ls -l
##cp -p dyn*.nc hist*.nc restart* gcm_out.txt $workdir/LMDZE
cp -p dyn*.nc hist*.nc restart* $workdir/LMDZE
