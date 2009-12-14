# This is a script in Bash.

# This script collects input files necessary for the program
# "etat0_lim" and runs "etat0_lim". Run this script in the directory
# where you want "etat0_lim" to run. The "ln" command does not report
# missing targets, so we check them first.

# Specify directories:

##in_dir=~/Documents/Around_LMDZ/Around_LMDZE/Input_etat0_lim
in_dir=$workdir/LMDZE/Input_etat0_lim

data_dir=~/Documents/Datasets

executable_dir=/usr/local/guez/LMDZ_work/LMDZE/Compil_prod_g95
##executable_dir=$workdir/LMDZE/Compil_prod

res_dir=/usr/local/guez/LMDZ_work/LMDZE/Results_etat0_lim
# (only for comparison with a previous run)

read -p \
    "Identifier for the set of parameters (\"*.def\" files and namelists)? " \
    iel_id

my_host=`hostname`

if [[ $my_host = brodie ]]
    then
    # (The signal "ERR" does not exist with the old Bash version
    # 2.05.8 on Brodie.)
    set -xe
else
    trap 'echo -e \\a; exit 1' ERR
    set -x
fi

if [[ $my_host = vierne ]]
    then
    # This is only useful if the program was compiled with G95:
    G95_FPU_UNDERFLOW=${G95_FPU_UNDERFLOW:+No}
fi

rm -f limit.nc start.nc startphy.nc coefoz_LMDZ.nc
# (If these are symbolic links then the Fortran program might not be
# able to replace them.)

test -d $in_dir/$iel_id

cd $data_dir
test -f Albedo.nc
test -f amipbc_sic_1x1.nc
test -f amipbc_sst_1x1.nc
##test -f ECPHY.nc
test -f ECDYN.nc
test -f landiceref.nc
test -f Relief.nc
test -f Ozone/Cariolle/Cariolle_v2_8/coefoz_v2_8.nc
test -f Rugos.nc
cd -

ln -s -f $data_dir/Albedo.nc $data_dir/amipbc_*.nc $data_dir/ECDYN.nc $data_dir/landiceref.nc $data_dir/Relief.nc $data_dir/Rugos.nc .

ln -sf $data_dir/Ozone/Cariolle/Cariolle_v2_8/coefoz_v2_8.nc coefoz.nc
ln -sf $data_dir/ECDYN.nc ECPHY.nc

rm -f *.def
test -f $in_dir/$iel_id/run.def
ln -s -f $in_dir/$iel_id/*.def .

if [[ $my_host = brodie ]]
    then
    rsh brodie01 \
    export F_PROGINF=YES \; \
    cd $PWD \; \
    $executable_dir/etat0_lim <$in_dir/$iel_id/etat0_lim_nml.txt \
    >etat0_lim_out.txt 2>etat0_lim_err.txt
elif [[ $my_host = zahir* ]]
    then
    hpmcount $executable_dir/etat0_lim <$in_dir/$iel_id/etat0_lim_nml.txt \
    >etat0_lim_out.txt 2>etat0_lim_err.txt
else
    time $executable_dir/etat0_lim <$in_dir/$iel_id/etat0_lim_nml.txt \
    >etat0_lim_out.txt 2>etat0_lim_err.txt
fi

echo -e '\a' # beep
set +x
read -p "Previous run number for comparison [none] ? "
if [[ -n $REPLY ]]
    then
    selective_diff.sh $res_dir/${REPLY# } .
fi
