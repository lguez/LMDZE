# This is a script in Bash.

trap 'exit 1' ERR
##set -x

read -p "Run number? " run_numb

read -p "Scalar or vector? (s/v) "
if [[ $REPLY = s ]]
    then
    read -p "Name of variable? " name
    read -p "Lower limit? " min
    read -p "Upper limit? " max
    read -p "Number of spatial dimensions? (2 or 3) " dim
    if ((dim == 2))
	then
	read -p "Average zonally? (y/n) "
	if [[ $REPLY = y ]]
	    then
            # Ferret script name: 
	    scr_name=rep_ave_pl.jnl
            # Parameters for the Ferret script:
	    set $name $run_numb $min $max
            # Name of output movie file:
	    out_name=${name}_x_ave.gif
	else
	    read -p "Interval between levels of color palette? " interv
            # Ferret script name: 
	    scr_name=rep_fill.jnl
            # Parameters for the Ferret script:
	    set $name $run_numb $min $max $interv
            # Name of output movie file:
	    out_name=$name.gif
	fi
    else
        # 3 spatial dimensions
	read -p "Interval between levels of color palette? " interv    
        # Ferret script name: 
	scr_name=rep_ave_f.jnl
        # Parameters for the Ferret script:
	set $name $run_numb $min $max $interv
        # Name of output movie file:
	out_name=${name}_x_ave.gif
    fi
else
    # Vector
    read -p "Name of variable 1? " name1
    read -p "Name of variable 2? " name2
    read -p "Standard vector length? " vec_len
    # Ferret script name: 
    scr_name=rep_vect.jnl
    # Parameters for the Ferret script:
    set $name1 $name2 $run_numb $vec_len
    # Name of output movie file:
    out_name=${name1}_${name2}.gif
fi

cd ~/Documents/Utilisation_LMDZ/Results_gcm/$run_numb
echo "Running \"$scr_name $*\""

movie.sh ~/Documents/Informatique_fonctionnement/Programs/LMDZE_program/Scripts/$scr_name "$*"
# (In older versions of Ferret, there should be a single argument
# after the script name. Therefore, for those versions, it is
# necessary to enclose $* between quotes.)

mv movie.gif $out_name
echo "Moved \"movie.gif\" to \"$out_name\"."
gifview -a $out_name &
